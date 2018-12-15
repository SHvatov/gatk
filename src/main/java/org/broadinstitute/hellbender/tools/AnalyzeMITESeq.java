package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.engine.ReferenceDataSource;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchMap;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Stream;

@DocumentedFeature
@CommandLineProgramProperties(
        summary = "(Experimental) Processes reads from a MITESeq experiment.",
        oneLineSummary = "(EXPERIMENTAL) Processes reads from a MITESeq experiment.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@BetaFeature
public class AnalyzeMITESeq extends GATKTool {
    @Argument(doc = "minimum quality score for analyzed portion of read", fullName = "min-q")
    private static int minQ = 30;

    @Argument(doc = "minimum size of analyzed portion of read", fullName = "min-length")
    private static int minLength = 15;

    @Argument(doc = "minimum number of wt calls flanking variant", fullName = "flanking-length")
    private static int flankingLength = 18;

    @Argument(doc = "reference indices of the ORF (1-based, closed), for example, '134-180,214-238'", fullName = "orf")
    private static String orfCoords;

    @Argument(doc = "minimum number of observations of reported variants", fullName = "min-variant-obs")
    private static long minVariantObservations = 8;

    @Argument(doc = "codon translation (a string of 64 amino acid codes", fullName = "codon-translation")
    private static String codonTranslation = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVZYZYSSSSZCWCLFLF";

    @Argument(doc = "output file prefix", fullName = "output-file-prefix", shortName = "O")
    private static String outputFilePrefix;

    @Argument(doc = "paired mode", fullName = "paired-mode")
    private static boolean pairedMode = true;

    private byte[] refSeq;
    private List<Interval> exonList;
    private List<Integer> refCodonValues;
    private final HopscotchMap<SNVCollectionCount, Long, SNVCollectionCount> variationCounts = new HopscotchMap<>(10000000);

    private long nReadsTotal = 0;
    private long nReadsUnmapped = 0;
    private long nReadsLowQuality = 0;
    private long nReadsWithLowQualityVariation = 0;
    private long[] refCoverage;
    private long[][] codonCounts;
    private GATKRead read1;

    private static final int LOWERCASE_MASK = 0xDF;
    private static final int N_REGULAR_CODONS = 64;
    private static final int FRAME_PRESERVING_INDEL_INDEX = 64;
    private static final int FRAME_SHIFTING_INDEL_INDEX = 65;
    private static final int CODON_COUNT_ROW_SIZE = 66;
    private static final int START_CODON = 0x0E;
    private static final int STOP_OCH = 0x30;
    private static final int STOP_AMB = 0x32;
    private static final int STOP_OPA = 0x38;
    private static final String[] labelForCodonValue = {
            "AAA","AAC","AAG","AAT","ACA","ACC","ACG","ACT","AGA","AGC","AGG","AGT","ATA","ATC","ATG","ATT",
            "CAA","CAC","CAG","CAT","CCA","CCC","CCG","CCT","CGA","CGC","CGG","CGT","CTA","CTC","CTG","CTT",
            "GAA","GAC","GAG","GAT","GCA","GCC","GCG","GCT","GGA","GGC","GGG","GGT","GTA","GTC","GTG","GTT",
            "TAA","TAC","TAG","TAT","TCA","TCC","TCG","TCT","TGA","TGC","TGG","TGT","TTA","TTC","TTG","TTT"
    };

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        initializeRefSeq();
        initializeExons();
        if ( codonTranslation.length() != N_REGULAR_CODONS ) {
            throw new UserException("codon-translation string must contain exactly 64 characters");
        }
    }

    @Override
    public void traverse() {
        // ignore unaligned reads and non-primary alignments
        final Stream<GATKRead> reads = getTransformedReadStream(ReadFilterLibrary.PRIMARY_LINE);

        if ( !pairedMode ) {
            reads.forEach(read -> {
                try {
                    final ReadReport readReport = processRead(read);
                    if ( readReport != null ) {
                        applyReport(readReport);
                    }
                } catch ( final Exception e ) {
                    throw new GATKException("Caught unexpected exception on read " +nReadsTotal + ": " + read.getName(), e);
                }
            });
        } else {
            reads.forEach(read -> {
                try {
                    if ( !read.isPaired() ) {
                        final ReadReport readReport = processRead(read);
                        if ( readReport != null ) {
                            applyReport(readReport);
                        }
                        return;
                    }
                    if ( read1 == null ) {
                        read1 = read;
                        return;
                    }
                    if ( !read1.getName().equals(read.getName()) ) {
                        System.out.println("Read " + read1.getName() + " has no mate.");
                        final ReadReport readReport = processRead(read1);
                        if ( readReport != null ) {
                            applyReport(readReport);
                        }
                        read1 = read;
                        return;
                    }
                    final ReadReport readReport1 = processRead(read1);
                    final ReadReport readReport2 = processRead(read);
                    if ( readReport1 == null ) {
                        if ( readReport2 != null ) {
                            applyReport(readReport2);
                        }
                    } else if ( readReport2 == null ) {
                        applyReport(readReport1);
                    } else {
                        applyReport(readReport1, readReport2);
                    }
                    read1 = null;
                } catch ( final Exception e ) {
                    throw new GATKException("Caught unexpected exception on read " +nReadsTotal + ": " + read.getName(), e);
                }
            });
        }
        if ( read1 != null ) {
            if ( read1.isPaired() ) {
                System.out.println("Read " + read1.getName() + " has no mate.");
            }
            try {
                final ReadReport readReport = processRead(read1);
                if ( readReport != null ) {
                    applyReport(readReport);
                }
            } catch ( final Exception e ) {
                throw new GATKException("Caught unexpected exception on read " +nReadsTotal + ": " + read1.getName(), e);
            }
        }
    }

    @Override
    public Object onTraversalSuccess() {
        writeVariationCounts(getVariationEntries());
        writeRefCoverage();
        writeCodonCounts();
        writeCodonFractions();
        writeAACounts();
        writeAAFractions();
        writeReadCounts();
        return null;
    }

    private List<SNVCollectionCount> getVariationEntries() {
        final long outputSize = variationCounts.stream().filter(entry -> entry.getCount() >= minVariantObservations).count();
        final List<SNVCollectionCount> variationEntries = new ArrayList<>((int)outputSize);
        for ( final SNVCollectionCount entry : variationCounts ) {
            if ( entry.getCount() > minVariantObservations ) {
                variationEntries.add(entry);
            }
        }
        variationEntries.sort((a,b) -> {
            int result = Long.compare(b.getCount(), a.getCount()); // descending order of count
            if ( result == 0 ) result = a.compareTo(b);
            return result;
        });
        return variationEntries;
    }

    private void writeVariationCounts( final List<SNVCollectionCount> variationEntries ) {
        final String variantsFile = outputFilePrefix + ".variantCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(variantsFile))) ) {
            for ( final SNVCollectionCount entry : variationEntries ) {
                writer.write(Long.toString(entry.getCount()));
                writer.write('\t');
                final List<SNV> snvs = entry.getSNVs();
                writer.write(Integer.toString(snvs.size()));
                String sep = "\t";
                for ( final SNV snv : snvs ) {
                    writer.write(sep);
                    sep = ", ";
                    writer.write(snv.toString());
                }
                writer.write(describeVariantsAsCodons(snvs));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+variantsFile, ioe);
        }
    }

    private String describeVariantsAsCodons( final List<SNV> snvs ) {
        final SortedMap<Integer, Integer> variantCodons = new TreeMap<>();
        for ( final SNV snv : snvs ) {

            //TODO: handle indels
            if ( snv.getRefCall() == '-' || snv.getVariantCall() == '-' ) continue;

            int codonIndex = 0;
            final Iterator<Interval> exonIterator = exonList.iterator();
            Interval currentInterval = exonIterator.next();
            while ( snv.getRefIndex() >= currentInterval.getEnd() ) {
                codonIndex += currentInterval.size();
                currentInterval = exonIterator.next();
            }
            if ( snv.getRefIndex() < currentInterval.getStart() ) {
                continue; // we're not in the currentInterval -- it's an intronic variation
            }
            codonIndex += snv.getRefIndex() - currentInterval.getStart();
            final int shift = 2 * (2 - (codonIndex % 3));
            codonIndex /= 3;
            int codonValue = variantCodons.getOrDefault(codonIndex, refCodonValues.get(codonIndex));
            final int callCode = "ACGT".indexOf(snv.getVariantCall());
            codonValue = (codonValue & ~(3 << shift)) | (callCode << shift);
            variantCodons.put(codonIndex, codonValue);
        }
        final StringBuilder sb = new StringBuilder();
        sb.append('\t').append(variantCodons.size());
        String prefix = "\t";
        for ( final Map.Entry<Integer, Integer> entry : variantCodons.entrySet() ) {
            final int codonValue = entry.getValue();
            sb.append(prefix).append(entry.getKey()).append(':')
                    .append(labelForCodonValue[refCodonValues.get(entry.getKey())])
                    .append('>').append(labelForCodonValue[codonValue]);
            prefix = ", ";
        }
        prefix = "\t";
        for ( final Map.Entry<Integer, Integer> entry : variantCodons.entrySet() ) {
            final char refAA = codonTranslation.charAt(refCodonValues.get(entry.getKey()));
            final int codonValue = entry.getValue();
            final boolean isStop = codonValue == STOP_OCH || codonValue == STOP_AMB || codonValue == STOP_OPA;
            final char altAA = codonTranslation.charAt(codonValue);
            sb.append(prefix).append(refAA==altAA ? 'S' : isStop ? 'N' : 'M').append(':')
                    .append(refAA).append('>').append(altAA);
            prefix = ", ";
        }
        return sb.toString();
    }

    private void writeRefCoverage() {
        final String refCoverageFile = outputFilePrefix + ".refCoverage";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(refCoverageFile))) ) {
            writer.write("RefPos\tCoverage\n");
            int refPos = 1;
            for ( final long coverage : refCoverage ) {
                writer.write(Integer.toString(refPos++));
                writer.write('\t');
                writer.write(Long.toString(coverage));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+refCoverageFile, ioe);
        }
    }

    private void writeCodonCounts() {
        final String codonCountsFile = outputFilePrefix + ".codonCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonCountsFile))) ) {
            final int nCodons = codonCounts.length;
            writer.write("AAA\tAAC\tAAG\tAAT\tACA\tACC\tACG\tACT\tAGA\tAGC\tAGG\tAGT\tATA\tATC\tATG\tATT\t"+
                    "CAA\tCAC\tCAG\tCAT\tCCA\tCCC\tCCG\tCCT\tCGA\tCGC\tCGG\tCGT\tCTA\tCTC\tCTG\tCTT\t"+
                    "GAA\tGAC\tGAG\tGAT\tGCA\tGCC\tGCG\tGCT\tGGA\tGGC\tGGG\tGGT\tGTA\tGTC\tGTG\tGTT\t"+
                    "TAA\tTAC\tTAG\tTAT\tTCA\tTCC\tTCG\tTCT\tTGA\tTGC\tTGG\tTGT\tTTA\tTTC\tTTG\tTTT\t"+
                    "NFS\tFS\tTotal\n");
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                long total = 0;
                for ( final long count : rowCounts ) {
                    writer.write(Long.toString(count));
                    writer.write('\t');
                    total += count;
                }
                writer.write(Long.toString(total));
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+codonCountsFile, ioe);
        }
    }

    private void writeCodonFractions() {
        final String codonFractFile = outputFilePrefix + ".codonFractions";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(codonFractFile))) ) {
            final int nCodons = codonCounts.length;
            writer.write("Codon  AAA  AAC  AAG  AAT  ACA  ACC  ACG  ACT  AGA  AGC  AGG  AGT  ATA  ATC  ATG  ATT"+
                    "  CAA  CAC  CAG  CAT  CCA  CCC  CCG  CCT  CGA  CGC  CGG  CGT  CTA  CTC  CTG  CTT"+
                    "  GAA  GAC  GAG  GAT  GCA  GCC  GCG  GCT  GGA  GGC  GGG  GGT  GTA  GTC  GTG  GTT"+
                    "  TAA  TAC  TAG  TAT  TCA  TCC  TCG  TCT  TGA  TGC  TGG  TGT  TTA  TTC  TTG  TTT"+
                    "  NFS   FS    Total\n");
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                writer.write(String.format("%5d", codonId + 1));
                final long[] rowCounts = codonCounts[codonId];
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : rowCounts ) {
                    writer.write(String.format("%5.1f", 100. * count / total));
                }
                writer.write(String.format(" %8d\n", total));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+codonFractFile, ioe);
        }
    }

    private void writeAACounts() {
        final String aaCountsFile = outputFilePrefix + ".aaCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(aaCountsFile))) ) {
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character,Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    String prefix = "";
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write(prefix);
                        prefix = "\t";
                        writer.write(chr);
                    }
                    writer.write('\n');
                }
                String prefix = "";
                for ( final long count : aaCounts.values() ) {
                    writer.write(prefix);
                    prefix = "\t";
                    writer.write(Long.toString(count));
                }
                writer.write('\n');
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+aaCountsFile, ioe);
        }
    }

    private void writeAAFractions() {
        final String aaFractFile = outputFilePrefix + ".aaFractions";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(aaFractFile))) ) {
            final int nCodons = codonCounts.length;
            for ( int codonId = 0; codonId != nCodons; ++codonId ) {
                final long[] rowCounts = codonCounts[codonId];
                final SortedMap<Character,Long> aaCounts = new TreeMap<>();
                for ( int codonValue = 0; codonValue != N_REGULAR_CODONS; ++codonValue ) {
                    aaCounts.merge(codonTranslation.charAt(codonValue), rowCounts[codonValue], Long::sum);
                }
                if ( codonId == 0 ) {
                    writer.write("Codon");
                    for ( final char chr : aaCounts.keySet() ) {
                        writer.write("    ");
                        writer.write(chr);
                    }
                    writer.write("    Total\n");
                }
                writer.write(String.format("%5d", codonId + 1));
                final long total = Arrays.stream(rowCounts).sum();
                for ( final long count : aaCounts.values() ) {
                    writer.write(String.format("%5.1f", 100. * count / total));
                }
                writer.write(String.format(" %8d\n", total));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+aaFractFile, ioe);
        }
    }

    private void writeReadCounts() {
        final String readCountsFile = outputFilePrefix + ".readCounts";
        try ( final OutputStreamWriter writer =
                      new OutputStreamWriter(new BufferedOutputStream(BucketUtils.createFile(readCountsFile))) ) {
            final DecimalFormat df = new DecimalFormat("0.000");
            writer.write("Total Reads:\t" + nReadsTotal + "\n");
            writer.write("Unmapped Reads:\t" + nReadsUnmapped + "\t" +
                    df.format(100.*nReadsUnmapped/nReadsTotal) + "%\n");
            writer.write("LowQ Reads:\t" + nReadsLowQuality + "\t" +
                    df.format(100.*nReadsLowQuality/nReadsTotal) + "%\n");
            writer.write("LowQ Var Reads:\t" + nReadsWithLowQualityVariation + "\t" +
                    df.format(100.*nReadsWithLowQualityVariation/nReadsTotal) + "%\n");
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write "+readCountsFile, ioe);
        }
    }

    // describes an interval as a pair of offsets (0-based, half-open) on the reference sequence.
    private final static class Interval {
        private final int start;
        private final int end;

        public Interval( final int begin, final int end ) {
            this.start = begin;
            this.end = end;
        }

        public int getStart() { return start; }
        public int getEnd() { return end; }
        public int size() { return end - start; }
    }

    private static final class SNV implements Comparable<SNV> {
        private final int refIndex;
        private final byte refCall;
        private final byte variantCall;

        public SNV( final int refIndex, final byte refCall, final byte variantCall ) {
            this.refIndex = refIndex;
            this.refCall = refCall;
            this.variantCall = variantCall;
        }

        public int getRefIndex() { return refIndex; }
        public byte getRefCall() { return refCall; }
        public byte getVariantCall() { return variantCall; }

        @Override
        public int hashCode() {
            return 47*(47*(47*refIndex + refCall) + variantCall);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNV && equals((SNV)obj);
        }

        public boolean equals( final SNV that ) {
            return this.refIndex == that.refIndex &&
                    this.refCall == that.refCall &&
                    this.variantCall == that.variantCall;
        }

        @Override
        public int compareTo( final SNV that ) {
            int result = Integer.compare(this.refIndex, that.refIndex);
            if ( result == 0 ) result = Byte.compare(this.refCall, that.refCall);
            if ( result == 0 ) result = Byte.compare(this.variantCall, that.variantCall);
            return result;
        }

        @Override
        public String toString() {
            return (refIndex+1) + ":" + (char)refCall + ">" + (char)variantCall;
        }
    }

    private static final class SNVCollectionCount
            implements Map.Entry<SNVCollectionCount, Long>, Comparable<SNVCollectionCount> {
        private static final SNV[] emptyArray = new SNV[0];
        private final SNV[] snvs;
        private long count;
        private final int hash;

        public SNVCollectionCount( final List<SNV> snvs ) {
            this.snvs = snvs.toArray(emptyArray);
            this.count = 1;
            int hashVal = 0;
            for ( final SNV snv : snvs ) {
                hashVal = 47*hashVal + snv.hashCode();
            }
            hash = 47*hashVal;
        }

        @Override
        public SNVCollectionCount getKey() { return this; }

        @Override
        public Long getValue() { return count; }

        @Override
        public Long setValue( final Long value ) {
            final Long result = count;
            count = value;
            return result;
        }

        public List<SNV> getSNVs() { return Arrays.asList(snvs); }
        public long getCount() { return count; }
        public void bumpCount() { count += 1; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof SNVCollectionCount && equals((SNVCollectionCount)obj);
        }

        public boolean equals( final SNVCollectionCount that ) {
            return this.hash == that.hash && Arrays.equals(this.snvs, that.snvs);
        }

        @Override
        public int hashCode() { return hash; }

        @Override
        public int compareTo( final SNVCollectionCount that ) {
            final int minSize = Math.min(this.snvs.length, that.snvs.length);
            int result = 0;
            for ( int idx = 0; idx != minSize; ++idx ) {
                result = this.snvs[idx].compareTo(that.snvs[idx]);
                if ( result != 0 ) break;
            }
            if ( result == 0 ) result = Integer.compare(this.snvs.length, that.snvs.length);
            return result;
        }
    }

    private static final class CodonTracker {
        private final Iterator<Interval> exonIterator;
        private Interval currentExon;
        private final int firstCodonIndex;
        private int codonPhase;
        private int codonValue;
        private int indelCount;
        private final List<Integer> codonValues;

        public CodonTracker( final List<Interval> exonList, final int refIndex ) {
            exonIterator = exonList.iterator();
            currentExon = exonIterator.next();
            int codonIndex = 0;
            while ( refIndex >= currentExon.getEnd() ) {
                codonIndex += currentExon.size();
                currentExon = exonIterator.next(); // there's a sentinel to keep us from going off the end
            }
            if ( refIndex > currentExon.getStart() ) {
                codonIndex += refIndex - currentExon.getStart();
            }
            codonPhase = codonIndex % 3;
            codonIndex /= 3;
            if ( codonPhase == 0 ) {
                firstCodonIndex = codonIndex;
            } else {
                firstCodonIndex = codonIndex + 1;
                codonPhase -= 3;
            }
            codonValue = 0;
            indelCount = 0;
            codonValues = new ArrayList<>();
        }

        public void push( final int refIndex, final byte call ) {
            if ( refIndex == currentExon.getEnd() ) {
                currentExon = exonIterator.next();
            }
            if ( refIndex < currentExon.getStart() ) {
                return;
            }

            final int callCode;
            switch ( call ) {
                case 'A': callCode = 0; break;
                case 'C': callCode = 1; break;
                case 'G': callCode = 2; break;
                case 'T': callCode = 3; break;
                case '+': indelCount += 1; return;
                case '-': indelCount -= 1; return;
                default: throw new GATKException("high quality call with value " + (char)call);
            }

            if ( indelCount != 0 ) {
                if ( codonPhase >= 0 ) {
                    codonValues.add(indelCount % 3 == 0 ? FRAME_PRESERVING_INDEL_INDEX : FRAME_SHIFTING_INDEL_INDEX);
                }

                // kill any further calling
                while ( exonIterator.hasNext() ) {
                    currentExon = exonIterator.next();
                }
                return;
            }

            codonValue = (codonValue << 2) | callCode;

            if ( ++codonPhase == 3 ) {
                final int val = codonValue & 0x3F;
                codonValues.add( val );
                codonPhase = 0;
                if ( val == STOP_OCH || val == STOP_AMB || val == STOP_OPA ) {
                    while ( exonIterator.hasNext() ) { // hit stop codon -- kill any further calling
                        currentExon = exonIterator.next();
                    }
                }
            }
        }

        public int getFirstCodonIndex() { return firstCodonIndex; }
        public List<Integer> getCodonValues() { return codonValues; }

        public boolean isConsistent( final CodonTracker that ) {
            final int startOverlap = Math.max(this.firstCodonIndex, that.firstCodonIndex);
            final int endOverlap = Math.min(this.getLastIndex(), that.getLastIndex());
            if ( startOverlap >= endOverlap ) return true;
            int idx1 = startOverlap - this.firstCodonIndex;
            int idx2 = startOverlap - that.firstCodonIndex;
            for ( int idx = startOverlap; idx != endOverlap; ++idx ) {
                if ( this.codonValues.get(idx1++).intValue() != that.codonValues.get(idx2++).intValue() ) {
                    return false;
                }
            }
            return true;
        }

        public void report( final long[][] codonCounts ) {
            int codonId = firstCodonIndex;
            for ( final int codonValue : codonValues ) {
                codonCounts[codonId++][codonValue] += 1;
            }
        }

        public void reportPair( final long[][] codonCounts, final CodonTracker that ) {
            report(codonCounts);
            if ( this.firstCodonIndex > that.firstCodonIndex ) {
                final int nnn = Math.min(that.codonValues.size(), this.firstCodonIndex - that.firstCodonIndex);
                int codonId = that.firstCodonIndex;
                for ( int idx = 0; idx != nnn; ++idx )
                    codonCounts[codonId++][that.codonValues.get(idx)] += 1;
            }
            if ( this.getLastIndex() < that.getLastIndex() ) {
                int nnn = Math.min(that.codonValues.size(), that.getLastIndex() - this.getLastIndex());
                int codonId = that.getLastIndex();
                int idx = that.codonValues.size();
                while ( nnn-- > 0 ) {
                    codonCounts[--codonId][that.codonValues.get(--idx)] += 1;
                }
            }
        }

        private int getLastIndex() { return firstCodonIndex + codonValues.size(); }
    }

    private final static class ReadReport {
        final List<Interval> refCoverage;
        final List<SNV> snvList;
        final CodonTracker codonTracker;

        public ReadReport( final List<Interval> refCoverage, final List<SNV> snvList, final CodonTracker codonTracker ) {
            this.refCoverage = refCoverage;
            this.snvList = snvList;
            this.codonTracker = codonTracker;
        }

        public List<Interval> getRefCoverage() { return refCoverage; }
        public List<SNV> getVariations() { return snvList; }
        public CodonTracker getCodonTracker() { return codonTracker; }
    }

    private void initializeRefSeq() {
        final ReferenceDataSource reference = ReferenceDataSource.of(referenceArguments.getReferencePath());
        final SAMSequenceDictionary seqDict = reference.getSequenceDictionary();
        if ( seqDict.size() != 1 ) {
            throw new UserException("Expecting a reference with a single contig. " +
                    "The supplied reference has " + seqDict.size() + " contigs.");
        }
        final SAMSequenceRecord tig0 = seqDict.getSequence(0);
        final int refSeqLen = tig0.getSequenceLength();
        final SimpleInterval wholeTig = new SimpleInterval(tig0.getSequenceName(), 1, refSeqLen);
        refSeq = Arrays.copyOf(reference.queryAndPrefetch(wholeTig).getBases(),refSeqLen);
        for ( int idx = 0; idx < refSeqLen; ++idx ) {
            switch ( refSeq[idx] &= LOWERCASE_MASK ) { // make into lower case
                case 'A': case 'C': case 'G': case 'T':
                    break;
                default:
                    throw new UserException("Reference sequence contains something other than A, C, G, and T.");
            }
        }
        refCoverage = new long[refSeq.length];
    }

    private void initializeExons() {
        exonList = new ArrayList<>();
        for ( final String coordPair : orfCoords.split(",") ) {
            final String[] coords = coordPair.split("-");
            if ( coords.length != 2 ) {
                throw new UserException("Can't interpret ORF as list of pairs of coords: " + orfCoords);
            }
            try {
                final int begin = Integer.valueOf(coords[0]);
                if ( begin < 1 ) {
                    throw new UserException("Coordinates of ORF are 1-based.");
                }
                final int end = Integer.valueOf(coords[1]);
                if ( end < begin ) {
                    throw new UserException("Found ORF end coordinate less than begin: " + orfCoords);
                }
                // convert 1-based, inclusive intervals to 0-based, half-open
                final Interval exon = new Interval(begin-1, end);
                exonList.add(exon);
            }
            catch ( final NumberFormatException nfe ) {
                throw new UserException("Can't interpret ORF coords as integers: " + orfCoords);
            }
            for ( int idx = 1; idx < exonList.size(); ++idx ) {
                if ( exonList.get(idx-1).getEnd() >= exonList.get(idx).getStart() ) {
                    throw new UserException("ORF coordinates are not sorted: " + orfCoords);
                }
            }
        }

        final int orfLen = exonList.stream().mapToInt(Interval::size).sum();
        if ( (orfLen % 3) != 0 ) {
            throw new UserException("ORF length must be divisible by 3.");
        }

        codonCounts = new long[orfLen / 3][];
        for ( int codonId = 0; codonId != codonCounts.length; ++codonId ) {
            codonCounts[codonId] = new long[CODON_COUNT_ROW_SIZE];
        }

        // it's helpful to have this 0-length sentinel at the end of the list
        exonList.add(new Interval(Integer.MAX_VALUE, Integer.MAX_VALUE));

        int refIndex = 0;
        final CodonTracker codonTracker = new CodonTracker(exonList, refIndex);
        final int refLen = refSeq.length;
        for ( int idx = 0; idx != refLen; ++idx ) {
            codonTracker.push(refIndex, refSeq[refIndex]);
            refIndex += 1;
        }
        refCodonValues = codonTracker.getCodonValues();
        if ( refCodonValues.size() != codonCounts.length ) {
            throw new UserException("Inconsistent parsing of reference sequence into expressed codons. " +
                                        "Do you have an upstream stop codon?");
        }
        if ( refCodonValues.get(0) != START_CODON ) {
            System.err.println("WARNING:  Your ORF does not start with ATG, as expected.");
        }
        final int lastCodon = refCodonValues.get(refCodonValues.size()-1);
        if ( lastCodon != STOP_OCH && lastCodon != STOP_AMB && lastCodon != STOP_OPA ) {
            System.err.println("WARNING:  Your ORF does not end with a stop codon, as expected.");
        }
    }

    private ReadReport processRead( final GATKRead read ) {
        nReadsTotal += 1;
        if (read.isUnmapped()) {
            nReadsUnmapped += 1;
            return null;
        }

        final byte[] quals = read.getBaseQualitiesNoCopy();

        // find initial end-trim
        int start = 0;
        int hiQCount = 0;
        while (start < quals.length) {
            if (quals[start] < minQ) {
                hiQCount = 0;
            } else if (++hiQCount == minLength) {
                break;
            }
            start += 1;
        }
        if (start == quals.length) {
            nReadsLowQuality += 1;
            return null;
        }
        start -= minLength - 1;

        // find final end-trim
        int end = quals.length - 1;
        hiQCount = 0;
        while (end >= 0) {
            if (quals[end] < minQ) {
                hiQCount = 0;
            } else if (++hiQCount == minLength) {
                break;
            }
            end -= 1;
        }
        end += minLength;

        return analyze(read, start, end);
    }

    private ReadReport analyze( final GATKRead read, final int start, final int end ) {
        final Cigar cigar = read.getCigar();

        // reads with a soft end-clip are no good unless the clip is "off the end" of the amplicon
        if ( cigar.getLastCigarElement().getOperator() == CigarOperator.S ) {
            if ( read.getEnd() != refSeq.length - 1 ) {
                nReadsWithLowQualityVariation += 1;
                return null;
            }
        }
        // reads with a soft start-clip are no good unless the clip is before the beginning of the amplicon
        if ( cigar.getFirstCigarElement().getOperator() == CigarOperator.S ) {
            if ( read.getStart() > 1 ) {
                nReadsWithLowQualityVariation += 1;
                return null;
            }
        }

        final Iterator<CigarElement> cigarIterator = cigar.getCigarElements().iterator();
        CigarElement cigarElement = cigarIterator.next();
        CigarOperator cigarOperator = cigarElement.getOperator();
        int cigarElementCount = cigarElement.getLength();

        final byte[] readSeq = read.getBasesNoCopy();
        final byte[] readQuals = read.getBaseQualitiesNoCopy();

        final List<SNV> variations = new ArrayList<>();

        int refIndex = read.getStart() - 1; // 0-based numbering
        int readIndex = 0;

        CodonTracker codonTracker = null;
        final List<Interval> refCoverageList = new ArrayList<>();
        int refCoverageBegin = 0;
        int refCoverageEnd = 0;

        while ( true ) {
            if ( readIndex >= start ) {
                if ( codonTracker == null ) {
                    codonTracker = new CodonTracker(exonList, refIndex);
                    refCoverageBegin = refIndex;
                    refCoverageEnd = refIndex;
                }
                if ( cigarOperator == CigarOperator.D ) {
                    variations.add(new SNV(refIndex, refSeq[refIndex], (byte)'-'));
                    codonTracker.push(refIndex, (byte)'-');
                } else if ( cigarOperator == CigarOperator.I ) {
                    // low-quality variations spoil the read
                    if ( readQuals[readIndex] < minQ ) {
                        nReadsWithLowQualityVariation += 1;
                        return null;
                    }
                    variations.add(new SNV(refIndex, (byte)'-', readSeq[readIndex]));
                    codonTracker.push(refIndex, (byte)'+');
                } else if ( cigarOperator == CigarOperator.M ) {
                    final byte call = (byte) (readSeq[readIndex] & LOWERCASE_MASK);
                    if (call != refSeq[refIndex]) {
                        // low-quality variations spoil the read
                        if (readQuals[readIndex] < minQ) {
                            nReadsWithLowQualityVariation += 1;
                            return null;
                        }
                        variations.add(new SNV(refIndex, refSeq[refIndex], readSeq[readIndex]));
                    }
                    if ( refIndex == refCoverageEnd ) {
                        refCoverageEnd += 1;
                    } else {
                        refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));
                        refCoverageBegin = refIndex;
                        refCoverageEnd = refIndex + 1;
                    }
                    codonTracker.push(refIndex, call);
                } else if ( cigarOperator != CigarOperator.S ) {
                    throw new GATKException("unanticipated cigar operator");
                }
            }

            if ( cigarOperator.consumesReadBases() ) {
                if ( ++readIndex == end ) {
                    break;
                }
            }

            if ( cigarOperator.consumesReferenceBases() ) {
                refIndex += 1;
            }

            if ( --cigarElementCount == 0 ) {
                cigarElement = cigarIterator.next();
                cigarOperator = cigarElement.getOperator();
                cigarElementCount = cigarElement.getLength();
            }
        }

        refCoverageList.add(new Interval(refCoverageBegin, refCoverageEnd));

        if ( !variations.isEmpty() ) {
            if ( refCoverageEnd - variations.get(variations.size() - 1).getRefIndex() < flankingLength ||
                    variations.get(0).getRefIndex() - refCoverageList.get(0).getStart() < flankingLength ) {
                variations.clear();
            }
        }
        return new ReadReport(refCoverageList, variations, codonTracker);
    }

    private void applyReport( final ReadReport readReport ) {

        readReport.getCodonTracker().report( codonCounts );

        final List<Interval> refCoverageList = readReport.getRefCoverage();
        for ( final Interval refInterval : refCoverageList ) {
            final int refIntervalEnd = refInterval.getEnd();
            for ( int idx = refInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                refCoverage[idx] += 1;
            }
        }

        reportVariations(readReport.getVariations());
    }

    private void reportVariations( final List<SNV> variations ) {
        if ( !variations.isEmpty() ) {
            final SNVCollectionCount newVal = new SNVCollectionCount(variations);
            final SNVCollectionCount oldVal = variationCounts.find(newVal);
            if ( oldVal != null ) {
                oldVal.bumpCount();
            } else {
                variationCounts.add(newVal);
            }
        }
    }

    private void applyReport( final ReadReport readReport1, final ReadReport readReport2 ) {
        if ( !readReport1.getCodonTracker().isConsistent(readReport2.getCodonTracker()) ) {
            return;
        }
        readReport1.getCodonTracker().reportPair(codonCounts, readReport2.getCodonTracker());

        final Iterator<Interval> refCoverageItr1 = readReport1.getRefCoverage().iterator();
        final Iterator<Interval> refCoverageItr2 = readReport2.getRefCoverage().iterator();
        Interval refInterval1 = refCoverageItr1.next();
        Interval refInterval2 = refCoverageItr2.next();
        Interval curInterval;
        if ( refInterval1.getStart() < refInterval2.getStart() ) {
            curInterval = refInterval1;
            refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
        } else {
            curInterval = refInterval2;
            refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
        }
        while ( refInterval1 != null || refInterval2 != null ) {
            final Interval testInterval;
            if ( refInterval1 == null ) {
                testInterval = refInterval2;
                refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
            } else if ( refInterval2 == null ) {
                testInterval = refInterval1;
                refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
            } else if ( refInterval1.getStart() < refInterval2.getStart() ) {
                testInterval = refInterval1;
                refInterval1 = refCoverageItr1.hasNext() ? refCoverageItr1.next() : null;
            } else {
                testInterval = refInterval2;
                refInterval2 = refCoverageItr2.hasNext() ? refCoverageItr2.next() : null;
            }
            if ( curInterval.getEnd() <= testInterval.getStart() ) {
                final int refIntervalEnd = curInterval.getEnd();
                for ( int idx = curInterval.getStart(); idx != refIntervalEnd; ++idx ) {
                    refCoverage[idx] += 1;
                }
                curInterval = testInterval;
            } else if ( testInterval.getEnd() > curInterval.getEnd() ) {
                curInterval = new Interval(curInterval.getStart(), testInterval.getEnd());
            }
        }
        final int refIntervalEnd = curInterval.getEnd();
        for ( int idx = curInterval.getStart(); idx != refIntervalEnd; ++idx ) {
            refCoverage[idx] += 1;
        }

        final List<SNV> variations1 = readReport1.getVariations();
        final List<SNV> variations2 = readReport2.getVariations();
        if ( variations1.isEmpty() ) {
            reportVariations(variations2);
        } else if ( variations2.isEmpty() ) {
            reportVariations(variations1);
        } else if ( new SNVCollectionCount(variations1).equals(new SNVCollectionCount(variations2)) ) {
            reportVariations(variations1);
        }
    }
}
