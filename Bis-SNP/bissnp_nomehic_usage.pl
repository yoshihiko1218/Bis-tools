#!/usr/bin/perl -w
## scNOMe-HiC / non-directional NOMe-seq wrapper for BisSNP.
##
## Companion to `bissnp_easy_usage.pl`, but tailored for non-directional
## bisulfite libraries (single-cell NOMe-HiC, scNMT-seq variants, etc.) where
## R2 reverse reads carry XG=CT and the BAM `is_reverse` flag does NOT
## reflect the bisulfite-conversion frame. Pairs with the XG-aware filter
## chain in BisSNP ≥1.1.
##
## Outputs (with default --prefix derived from BAM):
##   <prefix>.cyt.vcf                       cytosine VCF (raw, per-base methylation)
##   <prefix>.snp.vcf                       SNP VCF
##   <prefix>.GCH.txt                       per-read GCH methylation status
##   <prefix>.HCG.txt                       per-read HCG methylation status
##   <prefix>.cyt.filtered.sort.vcf         filtered + sorted cytosine VCF
##   <prefix>.cyt.filtered.sort.GCH.bed     per-base GCH methylation BED (per strand)
##   <prefix>.cyt.filtered.sort.HCG.bed     per-base HCG methylation BED (per strand)
##   (and GCG / HCH BEDs if --allC)

use strict;
use Getopt::Long;
use File::Basename;
use Cwd 'abs_path';

sub usage {
    print <<'USAGE';

Usage:
  perl bissnp_nomehic_usage.pl [Options] BISSNP_JAR INPUT_BAM REF_GENOME DBSNP_VCF

NON-DIRECTIONAL NOMe-seq / scNOMe-HiC wrapper. Runs BisulfiteGenotyper with
the flag set required for libraries where the XG tag (not BAM is_reverse)
defines the bisulfite-conversion frame, then post-processes the outputs
into per-base BED files. Per-read GCH/HCG methylation tables are emitted
inline by the genotyper.

PREREQUISITES
  - BisSNP jar version >= 1.1 (for XG-aware filter chain).
  - Java 8 on PATH or supplied via --java8.
  - INPUT_BAM coordinate-sorted, indexed, with MD tag. For scNOMe-HiC,
    BAM should already have Fix-D preprocessing applied (run
    scripts/strip_paired_all.py before this wrapper) — without it, the
    GATK 3.8 adaptor-clip filter silently drops the reverse mate of
    chimeric Hi-C ligations.
  - Helper scripts sortByRefAndCor.pl and vcf2bed6plus2.pl found in the
    same directory as this wrapper (default: same Bis-SNP/ install).

POSITIONAL ARGUMENTS
  BISSNP_JAR    Path to BisSNP jar (e.g. Bis-SNP.latest.jar or BisSNP-1.1.jar).
  INPUT_BAM     Input BAM. Must contain ReadGroup tag.
  REF_GENOME    Reference .fa with .fai sidecar.
  DBSNP_VCF     dbSNP VCF (e.g. Homo_sapiens_assembly38.dbsnp138.vcf).

OPTIONS
  --interval FILE   restrict analysis to a BED of intervals (e.g. chr22.bed)
                    or a single chromosome name passed via -L (e.g. chr22).
  --prefix STR      output filename prefix (default: BAM basename without .bam)
  --mem NUM         Java heap GB (default: 20)
  --nt NUM          BisSNP threads (default: 1; values > 1 are NOT safe with
                    NOMESEQ_MODE — keep at 1)
  --mmq NUM         min mapping quality (default: 30)
  --mbq NUM         min base quality (default: 5)
  --qual NUM        stand_call_conf (default: 1 — lowered from BisSNP default 20
                    because single-cell coverage is sparse)
  --mm NUM          max non-bisulfite mismatch fraction (default: 0.3)
  --minPatConv NUM  min CH-pattern conversion ratio (default: 1.0 = disabled)
  --minConv NUM     min 5'-converted Cs required (default: 0)
  --allC            also produce per-base BED for GCG and HCH
                    (default: only GCH and HCG)
  --skipBed         skip the BED post-processing chain (run BisSNP only)
  --java8 PATH      path to Java 8 binary
                    (default: /software/java/jdk1.8.0_191/bin/java)
  --dryRun          print commands without executing
  --help            show this message

EXAMPLE (Quest, scNOMe-HiC test cell)
  perl bissnp_nomehic_usage.pl \
      --interval chr22 --prefix bs_pr --mem 20 \
      /home/jmj7858/epifluidlab/software/Bis-tools/Bis-SNP/Bis-SNP.latest.jar \
      batch2.scD24.TAGCTT.dedup.calmd.allunpaired.bam \
      /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
      /gpfs/projects/b1198/epifluidlab/yoshii/reference/hg38/Homo_sapiens_assembly38.dbsnp138.vcf

USAGE
    exit 1;
}

my $interval   = "";
my $prefix     = "";
my $mem        = 20;
my $nt         = 1;
my $mmq        = 30;
my $mbq        = 5;
my $qual       = 1;
my $mm         = 0.3;
my $minPatConv = 1.0;
my $minConv    = 0;
my $allC       = 0;
my $skipBed    = 0;
my $java8      = "/software/java/jdk1.8.0_191/bin/java";
my $dryRun     = 0;
my $help       = 0;

GetOptions(
    "interval=s"   => \$interval,
    "prefix=s"     => \$prefix,
    "mem=i"        => \$mem,
    "nt=i"         => \$nt,
    "mmq=i"        => \$mmq,
    "mbq=i"        => \$mbq,
    "qual=i"       => \$qual,
    "mm=f"         => \$mm,
    "minPatConv=f" => \$minPatConv,
    "minConv=i"    => \$minConv,
    "allC"         => \$allC,
    "skipBed"      => \$skipBed,
    "java8=s"      => \$java8,
    "dryRun"       => \$dryRun,
    "help"         => \$help,
) or usage();

usage() if $help || scalar(@ARGV) != 4;

my ($BISSNP, $BAM, $REF, $DBSNP) = @ARGV;

foreach my $pair (["BISSNP_JAR",$BISSNP], ["INPUT_BAM",$BAM],
                  ["REF_GENOME",$REF], ["DBSNP_VCF",$DBSNP]) {
    my ($name, $f) = @$pair;
    die "$name not found: $f\n" unless -e $f;
}
die "REF .fai sidecar not found: ${REF}.fai\n" unless -e "${REF}.fai";

unless ($prefix) {
    $prefix = $BAM;
    $prefix =~ s|\.bam$||;
}

my $script_dir = dirname(abs_path(__FILE__));
# Helpers live in upstream Bis-tools/utils/; some forks copy them next to
# the wrapper. Search both locations.
sub find_helper {
    my $name = shift;
    foreach my $d ("$script_dir", "$script_dir/../utils") {
        my $p = "$d/$name";
        return $p if -e $p;
    }
    return undef;
}
my $SORT_VCF = find_helper("sortByRefAndCor.pl");
my $VCF2BED  = find_helper("vcf2bed6plus2.pl");

unless ($skipBed) {
    die "Helper sortByRefAndCor.pl not found near $script_dir\n" unless $SORT_VCF;
    die "Helper vcf2bed6plus2.pl not found near $script_dir\n"   unless $VCF2BED;
}

my $JAVA = "$java8 -Xmx${mem}G";

my $cyt_vcf      = "${prefix}.cyt.vcf";
my $snp_vcf      = "${prefix}.snp.vcf";
my $gch_per_read = "${prefix}.GCH.txt";
my $hcg_per_read = "${prefix}.HCG.txt";
my $cyt_sort     = "${prefix}.cyt.sort.vcf";
my $snp_sort     = "${prefix}.snp.sort.vcf";
my $cyt_filt     = "${prefix}.cyt.filtered.sort.vcf";

sub stamp { my $t = localtime(); return "[$t]"; }

sub run {
    my ($desc, $cmd) = @_;
    print STDERR stamp() . "  $desc\n$cmd\n";
    return if $dryRun;
    my $rc = system($cmd);
    die stamp() . "  FAILED ($desc): exit $rc\n" if $rc != 0;
}

sub step_genotyper {
    my $cmd = "$JAVA -jar $BISSNP \\\n"
            . "    -R $REF \\\n"
            . "    -I $BAM \\\n"
            . "    -D $DBSNP \\\n"
            . "    -T BisulfiteGenotyper \\\n"
            . "    -vfn1 $cyt_vcf \\\n"
            . "    -vfn2 $snp_vcf \\\n"
            . "    -gchreads $gch_per_read \\\n"
            . "    -cpgreads $hcg_per_read \\\n"
            . "    -out_modes NOMESEQ_MODE -sm GM \\\n"
            . "    -stand_call_conf $qual \\\n"
            . "    -nonDirectional -badMate \\\n"
            . "    -mm $mm -minPatConv $minPatConv -minConv $minConv \\\n"
            . "    -mmq $mmq -mbq $mbq -nt $nt";
    $cmd .= " \\\n    -L $interval" if $interval;
    run("BisulfiteGenotyper (NOMESEQ_MODE, non-directional)", $cmd);
}

sub step_sort {
    foreach my $vcf ($cyt_vcf, $snp_vcf) {
        my $out = $vcf;
        $out =~ s|\.vcf$|.sort.vcf|;
        run("Sort $vcf",
            "perl $SORT_VCF --k 1 --c 2 --tmp ./ $vcf ${REF}.fai > $out");
    }
}

sub step_filter {
    my $cmd = "$JAVA -jar $BISSNP \\\n"
            . "    -R $REF -T VCFpostprocess \\\n"
            . "    -C GCH -C HCH -C GCG -C HCG \\\n"
            . "    -oldVcf $cyt_sort -snpVcf $snp_sort \\\n"
            . "    -newVcf $cyt_filt \\\n"
            . "    -o $cyt_filt.cpgSummary.txt -minCT 1 -qual $qual";
    $cmd .= " \\\n    -L $interval" if $interval;
    run("VCFpostprocess (NOMe contexts)", $cmd);
}

sub step_to_bed {
    my @ctx = $allC ? qw(GCH HCG GCG HCH) : qw(GCH HCG);
    foreach my $c (@ctx) {
        run("vcf2bed6plus2 $c",
            "perl $VCF2BED --only_good_call --seperate_strand $cyt_filt $c");
    }
}

print STDERR stamp() . "  Output prefix: $prefix\n";

step_genotyper();

unless ($skipBed) {
    step_sort();
    step_filter();
    step_to_bed();
}

print STDERR stamp() . "  Done. Per-read tables: $gch_per_read, $hcg_per_read\n";
print STDERR stamp() . "  Per-base BEDs: ${prefix}.cyt.filtered.sort.{GCH,HCG}*.bed\n"
    unless $skipBed;
