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

public class FREECWorkflow extends OicrWorkflow {

    //Versions of tools
    private String freecVersion;
    private String samtoolsVersion;
    private String rModule;
    private String rLibDir;
    
    //FREEC
    private String chrLengthFile = "";
    private String templateType;
    private String freecVarCoeff = "";
    private String freecWindow;
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String  targetFile = "";
   
   
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
    private static final String FREEC_CV_DEFAULT         = "0.05";
    private static final String FREEC_WINDOW_DEFAULT_EX  = "500";
    private static final String FREEC_WINDOW_DEFAULT_WG  = "50000";
    private static final boolean DEFAULT_SKIP_IF_MISSING = true;  // Conditional provisioning
    private static final String FREEC_PREFIX = "freec_";
    private static final String RLIBDIR_BASE = "CNV.R_modules-";
    private final static String WG           = "WG";
    private final static String EX           = "EX";
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            this.normal        = getProperty("input_files_normal").split(",");
            this.tumor         = getProperty("input_files_tumor").split(",");
            this.chrLengthFile = getProperty("reference_len_file");
            this.templateType  = getProperty("template_type");

            
            String fileSkipFlag = this.getOptionalProperty("skip_missing_files", Boolean.toString(DEFAULT_SKIP_IF_MISSING));
            this.skipFlag = Boolean.valueOf(fileSkipFlag);
            
            //allow all template types other than WG (given that their have valid interval file associated)
            if (!this.templateType.equals(WG)) { 
                this.targetFile  = getProperty("target_file");
                this.freecWindow = getOptionalProperty("freec_window", FREEC_WINDOW_DEFAULT_EX);
            } else {
                this.freecWindow = getOptionalProperty("freec_window", FREEC_WINDOW_DEFAULT_WG);
            }
            
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length != 1) {
                        Logger.getLogger(FREECWorkflow.class.getName()).log(Level.SEVERE, "This workflow supports one normal sample only - "
                                + "check your .ini file");
                        return (null);
                    }
            }

            //=====================Application Versions           
            this.freecVersion    = getProperty("freec_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.rModule         = getProperty("R_module");
            this.rLibDir         = RLIBDIR_BASE + getProperty("Rlibs_version");
            this.freecVarCoeff   = getOptionalProperty("freec_var_coefficient", FREEC_CV_DEFAULT);
            
            //=============A special flag that determines if we need to sort/index
            String sortFlag = getProperty("do_sort");
            if (sortFlag == null || sortFlag.isEmpty()) {
                Logger.getLogger(FREECWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = sortFlag.isEmpty() || sortFlag.equalsIgnoreCase("false") ? false : true;
            }
            
            String manualFlag = getProperty("manual_output");
            if (manualFlag == null) {
                Logger.getLogger(FREECWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
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
            Logger.getLogger(FREECWorkflow.class.getName()).log(Level.SEVERE, null, ex);
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
            Logger.getLogger(FREECWorkflow.class.getName()).log(Level.WARNING, null, e);
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

          // We have this wf supporting one normal sample
          // but we use normal file index as id
          for (int n = 0; n < this.normal.length; n++) {
               for (int t = 0; t < this.tumor.length; t++) {
               if (this.templateType.equals(WG)) {
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
                             this.localInputTumorFiles[t], n + 1, null);
               } else if (!this.targetFile.isEmpty() && this.templateType.equals(EX)) {
                 // LAUNCH FREEC
                 launchFREEC(this.localInputNormalFiles[n],
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
     * FREEC configuring/launching
     */
    private void launchFREEC(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        // Job convertJob and create configFile
        Job freecJob = this.getWorkflow().createBashJob("freec_launch");
        String outputDir = this.dataDir + "FREEC." + id + "/";
        String resultID  = FREEC_PREFIX + this.makeBasename(inputTumor,  ".bam") + ".bam";
        freecJob.setCommand("module load " + this.rModule + ";"
                            + getWorkflowBaseDir() + "/dependencies/launchFREEC.pl"
                            + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --lenfile "      + this.chrLengthFile
                            + " --id "           + id
                            + " --freec "     + getWorkflowBaseDir() + "/bin/FREEC-" + this.freecVersion + "/freec"
                            + " --data-type " + this.templateType
                            + " --outdir "    + outputDir
                            + " --prefix "    + FREEC_PREFIX
                            + " --samtools "  + getWorkflowBaseDir() + "/bin/samtools-" + this.samtoolsVersion + "/samtools");
                            
        if (!this.freecVarCoeff.isEmpty()) {
         freecJob.getCommand().addArgument(" --var-coefficient " + this.freecVarCoeff);
        }
        if (!this.templateType.equals(WG)) {
         freecJob.getCommand().addArgument(" --target-file " + this.targetFile);
         freecJob.getCommand().addArgument(" --window "      + this.freecWindow);
        }
                
        freecJob.setMaxMemory("8000");
        if (parents != null) {
            for (Job p : parents) {
                freecJob.addParent(p);
            }
        }
        Log.stdout("Created FREEC launch Job");
        
        // Provision [tumor bam]_CNVs.p.value.txt, *_ratio.BedGraph, *_ratio_noNA.txt.png, *_sample.cpn, *_control.cpn
        SqwFile freecCNVFile = createOutputFile(outputDir + resultID + "_CNVs.p.value.txt", 
                                                "text/plain", this.manualOutput);
        freecCNVFile.setSkipIfMissing(skipFlag);
        freecCNVFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecCNVFile);
        
        SqwFile freecBedGraphFile = createOutputFile(outputDir + resultID + "_ratio.BedGraph",
                                                     "text/bed", this.manualOutput);
        freecBedGraphFile.setSkipIfMissing(skipFlag);
        freecBedGraphFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecBedGraphFile);
        
        SqwFile freecRatioPlotFile = createOutputFile(outputDir + resultID + "_ratio_noNA.txt.png",
                                                      "image/png", this.manualOutput);
        freecRatioPlotFile.setSkipIfMissing(skipFlag);
        freecRatioPlotFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecRatioPlotFile);
        
        // Raw (copy number profile) files
        SqwFile freecSampleCpnFile = createOutputFile(outputDir + resultID + "_sample.cpn",
                                                      "text/plain", this.manualOutput);
        freecSampleCpnFile.setSkipIfMissing(skipFlag);
        freecSampleCpnFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecSampleCpnFile);
        
        SqwFile freecControlCpnFile = createOutputFile(outputDir + this.makeBasename(inputNormal,  ".bam") + ".bam" 
                                                                           + "_control.cpn", "text/plain", this.manualOutput);
        freecControlCpnFile.setSkipIfMissing(skipFlag);
        freecControlCpnFile.getAnnotations().put("variation_calling_algorithm", "FREEC " + this.freecVersion);
        freecJob.addFile(freecControlCpnFile);
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
