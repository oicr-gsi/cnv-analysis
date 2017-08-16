package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

public class VarscanWorkflow extends OicrWorkflow {

    //Versions of tools
    private String varscanVersion;
    private String samtoolsVersion;
    private String rModule;
    private String rLibDir;
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String  refFasta;
    private String  templateType;

    
    /**   Varscan Filtering parameters:
     *    Defaults are set by the software, not this workflow
	--min-coverage	Minimum read depth at a position to make a call [8]
	--amp-threshold	Lower bound for log ratio to call amplification [0.25]
	--del-threshold	Upper bound for log ratio to call deletion (provide as positive number) [0.25]
	--min-region-size	Minimum size (in bases) for a region to be counted [10]
	--recenter-up	Recenter data around an adjusted baseline > 0 [0]
	--recenter-down	Recenter data around an adjusted baseline < 0 [0]
        
     */
    
    //Varscan
    private String varscanMinCoverage;
    private String varscanDelCoverage;
    private String varscanMinRegion;
    private String varscanRecenterUp;
    private String varscanRecenterDown;
    private String varscanPvalueThreshold;
    private String varscanJavaXmx;
    
    //Data
    private String[] normal;
    private String[] tumor;
    private String[] localInputNormalFiles;
    private String[] localInputTumorFiles;
    private String[] normalBases;
    private String[] tumorBases;

    //Misc
    private String dataDir;
    private String tempDir;

    private boolean skipFlag;
    private static final String RLIBDIR_BASE   = "CNV.R_modules-";
    private static final String VARSCAN_PREFIX = "varscan_";
    private final static String WG             = "WG";
    private final static String EX             = "EX";
    private final static String PVALUE         = "0.05";
    private static final boolean DEFAULT_SKIP_IF_MISSING = true;  // Conditional provisioning
    private static final String VARSCAN_JAVA_MEM = "4";
   
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            this.normal        = getProperty("input_files_normal").split(",");
            this.tumor         = getProperty("input_files_tumor").split(",");
            this.templateType  = getProperty("template_type");
           
            String fileSkipFlag = this.getOptionalProperty("skip_missing_files", Boolean.toString(DEFAULT_SKIP_IF_MISSING));
            this.varscanPvalueThreshold = this.getOptionalProperty("varscan_pvalue", PVALUE);
            this.varscanJavaXmx = this.getOptionalProperty("varscan_java_xmx", VARSCAN_JAVA_MEM);
            this.skipFlag = Boolean.valueOf(fileSkipFlag);           
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length != 1) {
                        Logger.getLogger(VarscanWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
            }

            //=====================Application Versions
            this.varscanVersion  = getProperty("varscan_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.refFasta        = getProperty("reference_fasta");
            this.rModule         = getProperty("R_module");
            this.rLibDir         = RLIBDIR_BASE + getProperty("Rlibs_version");
            this.varscanMinCoverage = getOptionalProperty("varscan_min_coverage", null);
            this.varscanDelCoverage = getOptionalProperty("varscan_del_coverage", null);
            this.varscanMinRegion   = getOptionalProperty("varscan_min_region", null);
            this.varscanRecenterUp  = getOptionalProperty("varscan_recenter_up", null);
            this.varscanRecenterUp  = getOptionalProperty("varscan_recenter_down", null);
            
            //=============A special flag that determines if we need to sort/index
            String sortFlag = getProperty("do_sort");
            if (sortFlag == null || sortFlag.isEmpty()) {
                Logger.getLogger(VarscanWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = sortFlag.isEmpty() || sortFlag.equalsIgnoreCase("false") ? false : true;
            }
            
            String manualFlag = getProperty("manual_output");
            if (manualFlag == null) {
                Logger.getLogger(VarscanWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
            } else {
                this.manualOutput = manualFlag.isEmpty() || manualFlag.equalsIgnoreCase("false") ? false : true;
            }

            // Register input files and set up local files
            this.localInputNormalFiles = new String[this.normal.length];
            this.localInputTumorFiles  = new String[this.tumor.length];
            this.normalBases           = new String[this.normal.length];
            this.tumorBases            = new String[this.tumor.length];
                   
            String[] types = {"normal", "tumor"};
            for (String type : types) {
                int listSize = type.equals("normal") ? this.normal.length : this.tumor.length;
                for (int fileIndex = 0; fileIndex < listSize; fileIndex++) {
                    String bamBasename;
                    /* vcftools analyzes filename, and if .gz is present, it will treat the file az bgzipped
                     * So, in our case we need to make sure that we don't have .gz in our files' names
                     */
                    if (type.equals("normal")) {
                        bamBasename = this.normal[fileIndex].substring(this.normal[fileIndex].lastIndexOf("/") + 1, this.normal[fileIndex].lastIndexOf(".bam"));
                        this.normalBases[fileIndex] = bamBasename.replaceAll(".gz.", ".");
                        this.normalBases[fileIndex] = this.normalBases[fileIndex].replaceAll(".fastq.annotated", "");
                    } else {
                        bamBasename = this.tumor[fileIndex].substring(this.tumor[fileIndex].lastIndexOf("/") + 1, this.tumor[fileIndex].lastIndexOf(".bam"));
                        this.tumorBases[fileIndex] = bamBasename.replaceAll(".gz.", ".");
                        this.tumorBases[fileIndex] = this.tumorBases[fileIndex].replaceAll(".fastq.annotated", "");
                    }

                    Log.stdout("CREATING FILE: input_bam_" + fileIndex + "_" + type);

                    SqwFile file = this.createFile("input_bam_" + fileIndex + "_" + type);
                    file.setType("application/bam");
                    file.setIsInput(true);
                    file.setForceCopy(false);

                    if (type.equals("normal")) {
                        file.setSourcePath(this.normal[fileIndex]);
                        this.localInputNormalFiles[fileIndex] = !this.doSort ? file.getProvisionedPath() : this.dataDir + bamBasename + ".sorted.bam";
                    } else {
                        file.setSourcePath(this.tumor[fileIndex]);
                        this.localInputTumorFiles[fileIndex] = !this.doSort ? file.getProvisionedPath() : this.dataDir + bamBasename + ".sorted.bam";
                    }
                }
            } // finished registering input bam files

            return this.getFiles();

        } catch (Exception ex) {
            Logger.getLogger(VarscanWorkflow.class.getName()).log(Level.SEVERE, null, ex);
            return (null);
        }
    }

    @Override
    public void setupDirectory() {
        try {
            this.dataDir = getProperty("data_dir");
            this.tempDir = "tempfiles/";
            if (this.dataDir == null || this.dataDir.isEmpty()) {
                this.dataDir = "data/";
            }
            if (!this.dataDir.endsWith("/")) {
                this.dataDir = this.dataDir.concat("/");
            }
            this.addDirectory(this.tempDir);
            this.addDirectory(this.dataDir);

        } catch (Exception e) {
            Logger.getLogger(VarscanWorkflow.class.getName()).log(Level.WARNING, null, e);
        }
    }

    @Override
    public void buildWorkflow() {

        try {
          
          List<Job> sortJobs = new ArrayList<Job>();

          if (this.doSort) {
          String[] types = {"normal", "tumor"};
            for (String type : types) {
                int listSize = type.equals("normal") ? this.normal.length : this.tumor.length;
                for (int fileIndex = 0; fileIndex < listSize; fileIndex++) {
                    String bamBasename = "";
                    String filePath = "";
                if (type.equals("normal")) {
                  bamBasename = this.makeBasename(this.localInputNormalFiles[fileIndex], ".bam");
                  filePath = this.normal[fileIndex];
                } else {
                  bamBasename = this.makeBasename(this.localInputTumorFiles[fileIndex], ".bam");
                  filePath = this.tumor[fileIndex];
                }
                  
                Job jobSamSort = this.getWorkflow().createBashJob("index_sort");
                jobSamSort.setCommand(getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools sort "
                                    + filePath + " "
                                    + this.dataDir + bamBasename); // localIndexed file
                jobSamSort.setMaxMemory("9000");
                 if (!this.queue.isEmpty()) {
                     jobSamSort.setQueue(this.queue);
                 }
                sortJobs.add(jobSamSort);
                }
           }
          }

           for (int n = 0; n < this.normal.length; n++) {
               for (int t = 0; t < this.tumor.length; t++) {

               /**
                * Need to think how to configure this for crosscheck properly 
                * if we don't have reference (TBD next iteration)
                */
               if (this.templateType.equals(WG) || this.templateType.equals(EX)) {
                 // LAUNCH BICseq
                 launchVarscan(this.localInputNormalFiles[n],
                              this.localInputTumorFiles[t], n + 1, sortJobs);
               
               } else {
                   throw new RuntimeException("Unsupported template type, workflow will terminate!");
               }
                   
              }
           }          
          // Summary job may be added in a next release

        } catch (Exception e) {
            Logger.getLogger(getClass().getName()).log(Level.SEVERE, null, e);
        }
    }
          

    /**
     * Varscan configuring/launching
     */
    private void launchVarscan(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        Job varscanJob   = this.getWorkflow().createBashJob("launch_varscan");   
        String outputDir = this.dataDir + "Varscan2." + id + "/"; 
        String resultID  = VARSCAN_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                         + this.makeBasename(inputTumor,  ".bam");
        varscanJob.setCommand("module load " + this.rModule + ";"
                            + getWorkflowBaseDir() + "/dependencies/launchVarscan2.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --output-dir "   + outputDir
                            + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                            + " --ref-fasta "    + refFasta
                            + " --java "         + getWorkflowBaseDir() + "/bin/jre" + getProperty("jre-version") + "/bin/java"
                            + " --varscan "      + getWorkflowBaseDir() + "/bin/VarScan.v" + varscanVersion + ".jar"
                            + " --id "           + resultID
                            + " --xmxmem "       + this.varscanJavaXmx
                            + " --p-value "      + this.varscanPvalueThreshold
                            + " --samtools "     + getWorkflowBaseDir() + "/bin/samtools-" 
                                                 + this.samtoolsVersion + "/samtools");
        if (null != this.varscanMinCoverage) {
            varscanJob.getCommand().addArgument(" --min-coverage " + this.varscanMinCoverage);
        }
        
        if (null != this.varscanDelCoverage) {
            varscanJob.getCommand().addArgument(" --del-coverage " + this.varscanDelCoverage);
        }
        
        if (null != this.varscanMinRegion) {
            varscanJob.getCommand().addArgument(" --min-region-size " + this.varscanMinRegion);
        }
                            
        if (null != this.varscanRecenterUp) {
            varscanJob.getCommand().addArgument(" --recenter-up " + this.varscanRecenterUp);
        }
        
        if (null != this.varscanRecenterDown) {
            varscanJob.getCommand().addArgument(" --recenter-down "+ this.varscanRecenterDown);
        }

        varscanJob.setMaxMemory("10000");
        if (parents != null) {
            for (Job p : parents) {
                varscanJob.addParent(p);
            }
        }
        Log.stdout("Created Varscan launch Job");
        
        // Provision files .copynumber, .copynumber.segmented, .copynumber.filtered
        //                 .copynumber.filtered.s_plot.png, .copynumber.filtered.s_plot.png       
        
        SqwFile varscanCopyFile = createOutputFile(outputDir + resultID + ".copynumber", 
                                                   "text/plain", this.manualOutput);
        varscanCopyFile.setSkipIfMissing(skipFlag);
        varscanCopyFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopyFile);
        
        SqwFile varscanCopySegFile = createOutputFile(outputDir + resultID + ".copynumber.segmented",
                                                      "text/plain", this.manualOutput);
        varscanCopySegFile.setSkipIfMissing(skipFlag);
        varscanCopySegFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopySegFile);
        
        SqwFile varscanCopyFilteredFile = createOutputFile(outputDir + resultID + ".copynumber.filtered",
                                                           "text/plain", this.manualOutput);
        varscanCopyFilteredFile.setSkipIfMissing(skipFlag);
        varscanCopyFilteredFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanCopyFilteredFile);
     
        SqwFile varscanWPlotFile = createOutputFile(outputDir + resultID + ".copynumber.filtered.w_plot.png",
                                                    "image/png", this.manualOutput);
        varscanWPlotFile.setSkipIfMissing(skipFlag);
        varscanWPlotFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanWPlotFile);
        
        SqwFile varscanSPlotFile = createOutputFile(outputDir + resultID + ".copynumber.filtered.s_plot.png",
                                                    "image/png", this.manualOutput);
        varscanSPlotFile.setSkipIfMissing(skipFlag);
        varscanSPlotFile.getAnnotations().put("variation_calling_algorithm", "Varscan " + this.varscanVersion);
        varscanJob.addFile(varscanSPlotFile);

    }
    
     
    /**
     * Utility function
     * 
     * @param path
     * @param extension
     * @return 
     */
    private String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }

}
