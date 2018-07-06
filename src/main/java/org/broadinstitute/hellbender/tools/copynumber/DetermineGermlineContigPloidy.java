package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.barclay.argparser.*;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineContigPloidyHybridADVIArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineContigPloidyModelArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.ContigCountDistributionCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.SimpleCount;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Determines the integer ploidy state of all contigs for germline samples given counts data. These should be either
 * HDF5 or TSV count files generated by {@link CollectReadCounts}.
 *
 * <h3>Introduction</h3>
 *
 * <p>Germline karyotyping is a frequently performed task in bioinformatics pipelines, e.g. for sex determination and
 * aneuploidy identification. This tool uses counts data for germline karyotyping. </p>
 *
 * <p>Performing germline karyotyping using counts data requires calibrating ("modeling") the technical coverage bias
 * and variance for each contig. The Bayesian model and the associated inference scheme implemented in
 * {@link DetermineGermlineContigPloidy} includes provisions for inferring and explaining away much of the technical
 * variation. Furthermore, karyotyping confidence is automatically adjusted for individual samples and contigs.</p>
 *
 * <p>Running {@link DetermineGermlineContigPloidy} is the first computational step in the GATK germline CNV calling
 * pipeline. It provides a baseline ("default") copy-number state for each contig/sample with respect to which the
 * probability of alternative states is allocated.</p>
 *
 * <h3>Python environment setup</h3>
 *
 * <p>The computation done by this tool, aside from input data parsing and validation, is performed outside of the Java
 * Virtual Machine and using the <em>gCNV computational python module</em>, namely {@code gcnvkernel}. It is crucial that
 * the user has properly set up a python conda environment with {@code gcnvkernel} and its dependencies
 * installed. If the user intends to run {@link DetermineGermlineContigPloidy} using one of the official GATK Docker images,
 * the python environment is already set up. Otherwise, the environment must be created and activated as described in the
 * main GATK README.md file.</p>
 *
 * <h3>Tool run modes</h3>
 *
 * <p>This tool has two operation modes as described below:</p>
 * <dl>
 *     <dt>COHORT mode:</dt>
 *
 *     <dd>If a ploidy model parameter path is not provided via the {@code model} argument, the tool will run in
 *     COHORT mode. In this mode, ploidy model parameters (e.g. coverage bias and variance for each contig) are
 *     inferred, along with baseline contig ploidy states of each sample. It is possible to run the tool over a subset
 *     of all intervals present in the input count files, which can be specified by -L; this can be used to pass a
 *     filtered interval list produced by {@link FilterIntervals} to mask intervals from modeling. The specified
 *     intervals must be present in all of the input count files.
 *
 *     <p>A TSV file specifying prior probabilities for each integer ploidy state and for each contig is required in this
 *     mode and must be specified via the {@code contig-ploidy-priors} argument. The following shows an example of
 *     such a table:</p>
 *     <br>
 *     <table border="1" width="80%">
 *         <tr>
 *             <td>CONTIG_NAME</td> <td>PLOIDY_PRIOR_0</td> <td>PLOIDY_PRIOR_1</td> <td>PLOIDY_PRIOR_2</td>
 *             <td>PLOIDY_PRIOR_3</td>
 *         </tr>
 *         <tr>
 *            <td>1</td> <td>0.01</td> <td>0.01</td> <td>0.97</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *            <td>2</td> <td>0.01</td> <td>0.01</td> <td>0.97</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *             <td>X</td> <td>0.01</td> <td>0.49</td> <td>0.49</td> <td>0.01</td>
 *         </tr>
 *         <tr>
 *            <td>Y</td> <td>0.50</td> <td>0.50</td> <td>0.00</td> <td>0.00</td>
 *         </tr>
 *     </table>
 *     <p>Note that the contig names appearing under {@code CONTIG_NAME} column must match contig names in the input
 *     counts files, and all contigs appearing in the input counts files must have a corresponding entry in the priors
 *     table. The order of contigs is immaterial in the priors table. The highest ploidy state is determined by the
 *     prior table (3 in the above example). A ploidy state can be strictly forbidden by setting its prior probability
 *     to 0. For example, the X contig in the above example can only assume 0 and 1 ploidy states.</p>
 *
 *     <p>The tool output in COHORT mode will contain two subdirectories, one ending with "-model" and the other
 *     ending with "-calls". The model subdirectory contains the inferred parameters of the ploidy model, which may
 *     be used later on for karyotyping one or more similarly-sequenced samples (see below).
 *
 *     The calls subdirectory contains one subdirectory for each sample, listing various sample-specific
 *     quantities such as the global read-depth, average ploidy, per-contig baseline ploidies, and per-contig
 *     coverage variance estimates.</p></dd>
 *
 *     <dt>CASE mode:</dt>
 *     <dd>If a path containing previously inferred ploidy model parameters is provided via the
 *     {@code model} argument, then the tool will run in CASE mode. In this mode, the parameters of the ploidy
 *     model are loaded from the provided directory and only sample-specific quantities are inferred. The modeled
 *     intervals are then specified by a file contained in the model directory, all interval-related arguments are
 *     ignored in this mode, and all model intervals must be present in all of the input count files. The tool output
 *     in CASE mode is only the "-calls" subdirectory and is organized similarly to that in COHORT mode.
 *
 *      <p>In CASE mode, the contig ploidy prior table is taken directly from the provided model parameters
 *      path and must be not provided again.</p></dd>
 * </dl>
 *
 * <h3>Important Remarks</h3>
 * <dl>
 *     <dt>Choice of hyperparameters:</dt>
 *     <dd><p>The quality of ploidy model parametrization and the sensitivity/precision of germline karyotyping are
 *     sensitive to the choice of model hyperparameters, including standard deviation of mean contig coverage bias
 *     (set using the {@code mean-bias-standard-deviation} argument), mapping error rate
 *     (set using the {@code mapping-error-rate} argument), and the typical scale of contig- and sample-specific
 *     unexplained variance (set using the {@code global-psi-scale} and {@code sample-psi-scale} arguments,
 *     respectively). It is crucial to note that these hyperparameters are <em>not</em> universal
 *     and must be tuned for each sequencing protocol and properly set at runtime.</p></dd>
 *
 *     <dt>Mosaicism and fractional ploidies:</dt>
 *     <dd><p>The model underlying this tool assumes integer ploidy states (in contrast to fractional/variable ploidy
 *     states). Therefore, it is to be used strictly on germline samples and for the purpose of sex determination,
 *     autosomal aneuploidy detection, or as a part of the GATK germline CNV calling pipeline. The presence of large somatic
 *     events and mosaicism (e.g., sex chromosome loss and somatic trisomy) will naturally lead to unreliable
 *     results. We strongly recommended inspecting genotyping qualities (GQ) from the tool output and considering to drop
 *     low-GQ contigs in downstream analyses. Finally, given the Bayesian status of this tool, we suggest including as many
 *     high-quality germline samples as possible for ploidy model parametrizaton in COHORT mode. This will downplay
 *     the role of questionable samples and will yield a more reliable estimation of genuine sequencing biases.</p></dd>
 *
 *     <dt>Coverage-based germline karyotyping:</dt>
 *     <dd><p>Accurate germline karyotyping requires incorporating SNP allele-fraction data and counts data in a
 *     unified probabilistic model and is beyond the scope of the present tool. The current implementation only uses
 *     counts data for karyotyping and while being fast, it may not provide the most reliable results.</p></dd>
 *
 * <h3>Usage examples</h3>
 *
 * <p>COHORT mode:</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy
 *   --run-mode COHORT \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --contig-ploidy-priors a_valid_ploidy_priors_table.tsv
 *   --output output_dir \
 *   --output-prefix normal_cohort
 * </pre>
 *
 * <p>COHORT mode (with optional interval filtering):</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy \
 *   --run-mode COHORT \
 *   -L intervals.interval_list \
 *   -XL blacklist_intervals.interval_list \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --contig-ploidy-priors a_valid_ploidy_priors_table.tsv
 *   --output output_dir \
 *   --output-prefix normal_cohort
 * </pre>
 *
 * <p>CASE mode:</p>
 * <pre>
 * gatk DetermineGermlineContigPloidy \
 *   --run-mode CASE \
 *   --model a_valid_ploidy_model_dir
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --output output_dir \
 *   --output-prefix normal_case
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Determines the baseline contig ploidy for germline samples given counts data",
        oneLineSummary = "Determines the baseline contig ploidy for germline samples given counts data",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class DetermineGermlineContigPloidy extends CommandLineProgram {
    public enum RunMode {
        COHORT, CASE
    }

    private static final String COHORT_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT = "cohort_determine_ploidy_and_depth.py";
    private static final String CASE_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT = "case_determine_ploidy_and_depth.py";

    //name of the interval file output by the python code in the model directory
    public static final String INPUT_MODEL_INTERVAL_FILE = "interval_list.tsv";

    public static final String MODEL_PATH_SUFFIX = "-model";
    public static final String CALLS_PATH_SUFFIX = "-calls";

    public static final String PLOIDY_STATE_PRIORS_FILE_LONG_NAME = "ploidy-state-priors";
    public static final String MAXIMUM_COUNT_LONG_NAME = "maximum-count";
    public static final String RUN_MODE_LONG_NAME = "run-mode";

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples.  " +
                    "Intervals must be identical and in the same order for all samples.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Tool run-mode.",
            fullName = RUN_MODE_LONG_NAME
    )
    private RunMode runMode;

    @Argument(
            doc = "Input file specifying ploidy-state priors.  This input is required in COHORT mode, " +
                    "but should not be provided in CASE mode.",
            fullName = PLOIDY_STATE_PRIORS_FILE_LONG_NAME,
            optional = true
    )
    private File inputPloidyStatePriorsFile;

    @Argument(
            doc = "Input ploidy-model directory.  This input is required in CASE mode, " +
                    "but should not be provided in COHORT mode.",
            fullName = CopyNumberStandardArgument.MODEL_LONG_NAME,
            optional = true
    )
    private String inputModelDir;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory for sample contig-ploidy calls and the contig-ploidy model parameters for " +
                    "future use.",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Argument(
            doc = "Maximum count to use in constructing coverage distributions.",
            fullName = MAXIMUM_COUNT_LONG_NAME,
            minValue = 0,
            optional = true
    )
    private int maximumCount = 250;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection
            = new OptionalIntervalArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineContigPloidyModelArgumentCollection germlineContigPloidyModelArgumentCollection =
            new GermlineContigPloidyModelArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineContigPloidyHybridADVIArgumentCollection germlineContigPloidyHybridADVIArgumentCollection =
            new GermlineContigPloidyHybridADVIArgumentCollection();

    private SimpleIntervalCollection specifiedIntervals;
    private File specifiedIntervalsFile;

    @Override
    protected void onStartup() {
        /* check for successful import of gcnvkernel */
        PythonScriptExecutor.checkPythonEnvironmentForPackage("gcnvkernel");
    }

    @Override
    protected Object doWork() {
        validateArguments();

        //read in count files and output intervals and contig x count distribution collections to temporary files
        final File contigCountDistributionCollectionsDir = IOUtils.createTempDir("contig-count-distribution-collections");
        final List<File> contigCountDistributionCollectionFiles =
                writeContigCountDistributionCollections(contigCountDistributionCollectionsDir, specifiedIntervals.getMetadata());

        //call python inference code
        final boolean pythonReturnCode = executeDeterminePloidyAndDepthPythonScript(
                contigCountDistributionCollectionFiles, specifiedIntervalsFile);

        if (!pythonReturnCode) {
            throw new UserException("Python return code was non-zero.");
        }

        logger.info("Germline contig ploidy determination complete.");

        return "SUCCESS";
    }

    private void validateArguments() {
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");

        if (inputModelDir != null) {
            Utils.validateArg(new File(inputModelDir).exists(),
                    String.format("Input ploidy-model directory %s does not exist.", inputModelDir));
            if (inputPloidyStatePriorsFile != null) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE mode, " +
                        "but ploidy-state priors were provided.");
            }
            if (intervalArgumentCollection.intervalsSpecified()) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE mode, " +
                        "but intervals were provided.");
            }
            //intervals are retrieved from the input model directory
            specifiedIntervalsFile = new File(inputModelDir, INPUT_MODEL_INTERVAL_FILE);
            IOUtils.canReadFile(specifiedIntervalsFile);
            specifiedIntervals = new SimpleIntervalCollection(specifiedIntervalsFile);
        } else {
            //get sequence dictionary and intervals from the first read-count file to use to validate remaining files
            //(this first file is read again below, which is slightly inefficient but is probably not worth the extra code)
            final File firstReadCountFile = inputReadCountFiles.get(0);
            final SimpleCountCollection firstReadCounts = SimpleCountCollection.read(firstReadCountFile);
            final SAMSequenceDictionary sequenceDictionary = firstReadCounts.getMetadata().getSequenceDictionary();
            final LocatableMetadata metadata = new SimpleLocatableMetadata(sequenceDictionary);

            if (intervalArgumentCollection.intervalsSpecified()) {
                logger.info("Intervals specified...");
                CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);
                specifiedIntervals = new SimpleIntervalCollection(metadata,
                        intervalArgumentCollection.getIntervals(sequenceDictionary));
            } else {
                logger.info(String.format("Retrieving intervals from first read-count file (%s)...",
                        firstReadCountFile));
                specifiedIntervals = new SimpleIntervalCollection(metadata, firstReadCounts.getIntervals());
            }

            //in cohort mode, intervals are specified via -L; we write them to a temporary file
            specifiedIntervalsFile = IOUtils.createTempFile("intervals", ".tsv");
            specifiedIntervals.write(specifiedIntervalsFile);
        }

        if (runMode.equals(RunMode.COHORT)) {
            logger.info("Running the tool in COHORT mode...");
            Utils.validateArg(inputReadCountFiles.size() > 1, "At least two samples must be provided in " +
                    "COHORT mode");
            if (inputModelDir != null) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in COHORT mode, " +
                        "but ploidy-model directory was provided.");
            }
            if (inputPloidyStatePriorsFile == null) {
                throw new UserException.BadInput("Ploidy-state priors must be provided in COHORT mode.");
            }
            IOUtils.canReadFile(inputPloidyStatePriorsFile);
        } else { // case run-mode
            logger.info("Running the tool in CASE mode...");
            Utils.validateArg(inputModelDir != null, "An input ploidy-model directory must be provided in " +
                    "CASE mode.");
            Utils.validateArg(new File(inputModelDir).exists(),
                    String.format("Input ploidy-model directory %s does not exist.", inputModelDir));
            if (inputPloidyStatePriorsFile != null) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE mode, " +
                        "but ploidy-state priors were provided.");
            }
            if (intervalArgumentCollection.intervalsSpecified()) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE mode, " +
                        "but intervals were provided.");
            }

            //get sequence dictionary and intervals from the first read-count file to use to validate remaining files
            //(this first file is read again below, which is slightly inefficient but is probably not worth the extra code)
            final File firstReadCountFile = inputReadCountFiles.get(0);
            specifiedIntervals = CopyNumberArgumentValidationUtils.resolveIntervals(
                    firstReadCountFile, intervalArgumentCollection, logger);
        }

        Utils.nonNull(outputPrefix);
        ParamUtils.isPositiveOrZero(maximumCount, "Maximum count must be non-negative.");
        Utils.validateArg(new File(outputDir).exists(),
                String.format("Output directory %s does not exist.", outputDir));

        germlineContigPloidyModelArgumentCollection.validate();
        germlineContigPloidyHybridADVIArgumentCollection.validate();
    }

    private List<File> writeContigCountDistributionCollections(final File contigCountDistributionCollectionsDir,
                                                               final LocatableMetadata metadata) {
        logger.info("Validating and aggregating per-contig count distributions from input read-count files...");
        final int numSamples = inputReadCountFiles.size();
        final List<File> contigCountDistributionCollectionFiles = new ArrayList<>(numSamples);
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        final Set<SimpleInterval> intervalSubset = new HashSet<>(specifiedIntervals.getRecords());

        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-count file %s (%d / %d)",
                    inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            if (!CopyNumberArgumentValidationUtils.isSameDictionary(
                    readCounts.getMetadata().getSequenceDictionary(), metadata.getSequenceDictionary())) {
                logger.warn("Sequence dictionary for read-count file %s does not match that " +
                        "in other read-count files.", inputReadCountFile);
            }
            final List<SimpleCount> subsetReadCounts = readCounts.getRecords().stream()
                    .filter(c -> intervalSubset.contains(c.getInterval()))
                    .collect(Collectors.toList());
            Utils.validateArg(subsetReadCounts.size() == intervalSubset.size(),
                    String.format("Intervals for read-count file %s do not contain all specified intervals.",
                            inputReadCountFile));
            //calculate per-contig count distributions and write temporary file for this sample
            final File contigCountDistributionCollectionFile =
                    new File(contigCountDistributionCollectionsDir, String.format("SAMPLE-%d.tsv", sampleIndex));
            new ContigCountDistributionCollection(readCounts, intervalSubset, maximumCount)
                    .write(contigCountDistributionCollectionFile);
            contigCountDistributionCollectionFiles.add(contigCountDistributionCollectionFile);
        }
        return contigCountDistributionCollectionFiles;
    }

    private boolean executeDeterminePloidyAndDepthPythonScript(final List<File> contigCountDistributionCollectionFiles,
                                                               final File intervalsFile) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final String outputDirArg = Utils.nonEmpty(outputDir).endsWith(File.separator)
                ? outputDir
                : outputDir + File.separator;    //add trailing slash if necessary
        final List<String> arguments = new ArrayList<>(Collections.singletonList(
                "--output_calls_path=" + outputDirArg + outputPrefix + CALLS_PATH_SUFFIX));
        arguments.add("--contig_count_distribution_collection_files");
        arguments.addAll(contigCountDistributionCollectionFiles.stream().map(File::getAbsolutePath).collect(Collectors.toList()));
        arguments.addAll(germlineContigPloidyModelArgumentCollection.generatePythonArguments(runMode));
        arguments.addAll(germlineContigPloidyHybridADVIArgumentCollection.generatePythonArguments());

        final String script;
        if (runMode == RunMode.COHORT) {
            script = COHORT_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT;
            arguments.add("--interval_list=" + intervalsFile.getAbsolutePath());
            arguments.add("--ploidy_state_priors_table=" + inputPloidyStatePriorsFile.getAbsolutePath());
            arguments.add("--output_model_path=" + outputDirArg + outputPrefix + MODEL_PATH_SUFFIX);
        } else {
            script = CASE_DETERMINE_PLOIDY_AND_DEPTH_PYTHON_SCRIPT;
            arguments.add("--input_model_path=" + inputModelDir);
        }
        return executor.executeScript(
                new Resource(script, GermlineCNVCaller.class),
                null,
                arguments);
    }
}
