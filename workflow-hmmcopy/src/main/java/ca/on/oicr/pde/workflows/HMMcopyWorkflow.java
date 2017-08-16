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

public class HMMcopyWorkflow extends OicrWorkflow {

    //Versions of tools
    private String samtoolsVersion;
    private String hmmcopyVersion;
    private String rModule;
    private String rLibDir;
    
    private String templateType = "";
            
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String[] supportedChromosomes;
    
    //HMMcopy
    private String refGCfile;
    private String refMAPfile;
  
    
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
    private static final boolean DEFAULT_SKIP_IF_MISSING = true;  // Conditional provisioning
    private static final String HMMCOPY_PREFIX = "hmmcopy_";
    private static final String RLIBDIR_BASE   = "CNV.R_modules-";
    private final static String WG             = "WG";
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            this.normal        = getProperty("input_files_normal").split(",");
            this.tumor         = getProperty("input_files_tumor").split(",");
            this.refGCfile     = getProperty("reference_gc");
            this.refMAPfile    = getProperty("reference_map");
            this.templateType  = getProperty("template_type");
            this.supportedChromosomes = getProperty("supported_chromosomes").split(",");

            
            String fileSkipFlag = this.getOptionalProperty("skip_missing_files", Boolean.toString(DEFAULT_SKIP_IF_MISSING));
            this.skipFlag = Boolean.valueOf(fileSkipFlag);
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length != 1) {
                        Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
            }
            
            //=====================Issue a warning if have anything other than WG, the wf will terminate shortly after
            if (this.templateType == null && !this.templateType.equals(WG)) {
                 Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.SEVERE, "template type set to {0} which is not supported for this workflow", this.templateType);
            }

            //=====================Application Versions           
            this.samtoolsVersion = getProperty("samtools_version");
            this.hmmcopyVersion  = getProperty("hmmcopy_version");
            this.rModule         = getProperty("R_module");
            this.rLibDir         = RLIBDIR_BASE + getProperty("Rlibs_version");
            
            //=============A special flag that determines if we need to sort/index
            String sortFlag = getProperty("do_sort");
            if (sortFlag == null || sortFlag.isEmpty()) {
                Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = sortFlag.isEmpty() || sortFlag.equalsIgnoreCase("false") ? false : true;
            }
            
            String manualFlag = getProperty("manual_output");
            if (manualFlag == null) {
                Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
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
            Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.SEVERE, null, ex);
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
            Logger.getLogger(HMMcopyWorkflow.class.getName()).log(Level.WARNING, null, e);
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
                * TODO need to think how to configure this for crosscheck properly 
                * if we don't have reference (TBD next iteration)
                */
               if (this.templateType.equals(WG)) {
                 // LAUNCH HMMcopy
                 launchHMMcopy(this.localInputNormalFiles[n],
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
     * HMMcopy configuring/launching
     */
    private void launchHMMcopy(String inputNormal, String inputTumor, int id, List<Job> parents) {
        
        String[] allInputs = {inputNormal, inputTumor};
        String outputDir = this.dataDir + "HMMcopy." + id + "/";
        List<Job> setupJobs = new ArrayList<Job>();
       
        // Inputs converted into .wig and then - .bw format
        for (String inFile : allInputs) {
        // =============== indexing =========================================
        Job indexJob = this.getWorkflow().createBashJob("hmmcopy_index");
        indexJob.setCommand("mkdir -p " + outputDir + ";"
                          + getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter -b " + inFile);
        indexJob.setMaxMemory("4000");
            
        if (parents != null) {
          for (Job p : parents) {
             indexJob.addParent(p);
          }
        }
            
        //============ converting to wig format =============================
        Job convertJob = this.getWorkflow().createBashJob("hmmcopy_convert");
        convertJob.setCommand(getWorkflowBaseDir() + "/dependencies/convertHMMcopy.pl "
                            + " --read-counter " + getWorkflowBaseDir() + "/bin/HMMcopy-" + this.hmmcopyVersion + "/bin/readCounter "
                            + " --input "  + inFile
                            + " --output " + this.makeBasename(inFile, ".bam") + "_reads.wig");
        convertJob.setMaxMemory("4000");
        convertJob.addParent(indexJob);
            
        setupJobs.add(convertJob);
        }
       
       // Launch HMMcopy scripts, provision results
       //================== run HMMcopy ====================================
        Job hmmJob = this.getWorkflow().createBashJob("hmmcopy_launch");
        String resultID = HMMCOPY_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                         + this.makeBasename(inputTumor,  ".bam");
        hmmJob.setCommand("module load " + this.rModule + ";"
                        + getWorkflowBaseDir() + "/dependencies/launchHMMcopy.pl "
                        + " --r-libdir "     + getWorkflowBaseDir() + "/bin/" + this.rLibDir
                        + " --normal-wig "   + this.makeBasename(inputNormal, ".bam") + "_reads.wig "
                        + " --tumor-wig "    + this.makeBasename(inputTumor, ".bam") + "_reads.wig "
                        + " --cg-file "      + this.refGCfile
                        + " --map-file "     + this.refMAPfile
                        + " --hmm-script "   + getWorkflowBaseDir() + "/dependencies/run_HMMcopy.r"
                        + " --output-base "  + outputDir + resultID);

        hmmJob.setMaxMemory("6000");
        for (Job p : setupJobs) {
            hmmJob.addParent(p);
        }
        
        Log.stdout("Created HMMcopy launch Job");
        
        // Provision .seg, .tsv, .bias_plot.png, .c_plot.chr*.png, .s_plot.chr*.png
        SqwFile hmmcopySegFile = createOutputFile(outputDir + resultID + ".seg", "text/plain", this.manualOutput);
        hmmcopySegFile.setSkipIfMissing(skipFlag);
        hmmcopySegFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopySegFile);
        
        SqwFile hmmcopyTsvFile = createOutputFile(outputDir + resultID + ".tsv", "text/plain", this.manualOutput);
        hmmcopyTsvFile.setSkipIfMissing(skipFlag);
        hmmcopyTsvFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopyTsvFile);
        
        SqwFile hmmcopyBiasPlotFile = createOutputFile(outputDir + resultID + ".png", "image/png", this.manualOutput);
        hmmcopyBiasPlotFile.setSkipIfMissing(skipFlag);
        hmmcopyBiasPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
        hmmJob.addFile(hmmcopyBiasPlotFile);
        
        for(String chrom : this.supportedChromosomes) {

            SqwFile hmmcopyCPlotFile = createOutputFile(outputDir + resultID + ".c_plot." + chrom + ".png", "image/png", this.manualOutput);
            hmmcopyCPlotFile.setSkipIfMissing(skipFlag);
            hmmcopyCPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
            hmmJob.addFile(hmmcopyCPlotFile);
            
            SqwFile hmmcopySPlotFile = createOutputFile(outputDir + resultID + ".s_plot." + chrom + ".png", "image/png", this.manualOutput);
            hmmcopySPlotFile.setSkipIfMissing(skipFlag);
            hmmcopySPlotFile.getAnnotations().put("variation_calling_algorithm", "HMMcopy " + this.hmmcopyVersion);
            hmmJob.addFile(hmmcopySPlotFile);
   
        }
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
