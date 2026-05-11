/**
 * XG-aware BisulfiteMismatchReadsFilter (full fix).
 *
 * Consults the XG tag (set by bhmem-XG) to decide which bisulfite conversion
 * pattern is tolerated for this read:
 *   XG=CT  → tolerate  ref=C  read=T   (read came from + strand template)
 *   XG=GA  → tolerate  ref=G  read=A   (read came from - strand template)
 * Falls back to the 2018 simpleReverseComplement + library-agnostic logic
 * only if the XG tag is missing.
 */
package edu.usc.epigenome.uecgatk.bissnp.filters;

import edu.usc.epigenome.uecgatk.bissnp.BaseUtilsMore;
import edu.usc.epigenome.uecgatk.bissnp.BisSNPUtils;
import edu.usc.epigenome.uecgatk.bissnp.BisulfiteSAMConstants;
import htsjdk.samtools.SAMRecord;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.engine.filters.ReadFilter;
import org.broadinstitute.gatk.utils.BaseUtils;

public class BisulfiteMismatchReadsFilter extends ReadFilter {

	@Argument(fullName = "max_mismatches", shortName = "mm",
			doc = "Maximum percentage of non-bisulfite mismatches within a read for a read to be used for calling. Default: 0.3",
			required = false)
	public static double MAX_MISMATCHES = 0.3;

	@Override
	public boolean filterOut(SAMRecord read) {
		if (read.getStringAttribute(BisulfiteSAMConstants.MD_TAG) == null) return false;
		try {
			return hasTooManyBisulfiteMismathces(read);
		} catch (Exception e) { e.printStackTrace(); }
		return false;
	}

	public static boolean hasTooManyBisulfiteMismathces(SAMRecord read) throws Exception {
		byte[] refBases = BaseUtilsMore.toUpperCase(
				BisSNPUtils.modifyRefSeqByCigar(BisSNPUtils.refStrFromMd(read), read.getCigarString()));
		byte[] bases    = BaseUtilsMore.toUpperCase(BisSNPUtils.getClippedReadsBase(read));

		String xg = read.getStringAttribute("XG");
		byte tolRef = 0, tolBase = 0;
		if ("CT".equals(xg)) { tolRef = (byte)'C'; tolBase = (byte)'T'; }
		else if ("GA".equals(xg)) { tolRef = (byte)'G'; tolBase = (byte)'A'; }

		if (tolRef == 0) {
			// XG missing: fall back to 2018 simpleReverseComplement + library-agnostic
			boolean neg = read.getReadNegativeStrandFlag();
			if (neg) {
				bases = BaseUtils.simpleReverseComplement(bases);
				refBases = BaseUtils.simpleReverseComplement(refBases);
			}
			int len = Math.min(bases.length, refBases.length);
			int nmm = 0;
			for (int i = 0; i < len; i++) {
				if (BaseUtils.basesAreEqual(refBases[i], bases[i])) continue;
				if (BaseUtilsMore.isBisulfiteMismatch(refBases[i], bases[i], neg)) nmm++;
				if (nmm > MAX_MISMATCHES * len) return true;
			}
			return false;
		}

		// XG present: compare in the read's aligned frame directly; tolerate the
		// one conversion pattern indicated by XG.
		int len = Math.min(bases.length, refBases.length);
		int nmm = 0;
		for (int i = 0; i < len; i++) {
			if (BaseUtils.basesAreEqual(refBases[i], bases[i])) continue;
			if (refBases[i] == tolRef && bases[i] == tolBase) continue;
			nmm++;
			if (nmm > MAX_MISMATCHES * len) return true;
		}
		return false;
	}

	@Override
	public boolean filterOut(SAMRecord a, SAMRecord b) { return false; }
}
