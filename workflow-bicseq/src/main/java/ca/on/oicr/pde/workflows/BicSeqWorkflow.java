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

public class BicSeqWorkflow extends OicrWorkflow {

    //Versions of tools
    private String bicseqVersion;
    private String samtoolsVersion;
    private String rModule;
               
    //References
    private boolean manualOutput;
    private boolean doSort = true;
    private String  queue;
    private String  templateType = "";
       
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
    private int bicseqInterval;
    private int bicseqSpread;

    private boolean skipFlag;
    private static final String BICSEQ_I_DEFAULT         = "150";
    private static final String BICSEQ_S_DEFAULT         = "20";
    private static final boolean DEFAULT_SKIP_IF_MISSING = true;  // Conditional provisioning
    private static final String BICSEQ_PREFIX  = "bicseq_";
    private final static String WG           = "WG";    
    
    @Override
    public Map<String, SqwFile> setupFiles() {

        try {

            this.normal        = getProperty("input_files_normal").split(",");
            this.tumor         = getProperty("input_files_tumor").split(",");
            this.templateType  = getProperty("template_type");
            
            String fileSkipFlag = this.getOptionalProperty("skip_missing_files", Boolean.toString(DEFAULT_SKIP_IF_MISSING));
            this.skipFlag = Boolean.valueOf(fileSkipFlag);
                       
            this.queue = getProperty("queue");
            
            if (this.tumor.length != this.normal.length) {
                    if (this.normal.length != 1) {
                        Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.SEVERE, "numbers for normal and tumor bam files are not the same, crosscheck isn't forced - "
                                + "check your .ini file");
                        return (null);
                    }
            }

            //=====================Issue a warning if have anything other than WG, the wf will terminate shortly after
            if (this.templateType == null && !this.templateType.equals(WG)) {
                 Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.SEVERE, "template type set to {0} which is not supported for this workflow", this.templateType);
            }
            
            //=====================Application Versions           
            this.bicseqVersion   = getProperty("bicseq_version");
            this.samtoolsVersion = getProperty("samtools_version");
            this.rModule         = getProperty("R_module");
            
            //=============A special flag that determines if we need to sort/index
            String sortFlag = getProperty("do_sort");
            if (sortFlag == null || sortFlag.isEmpty()) {
                Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.WARNING, "do_sort is not set, will deduce it from the names of input files");
                this.doSort = !this.normal[0].contains("sorted");
            } else {
                this.doSort = sortFlag.isEmpty() || sortFlag.equalsIgnoreCase("false") ? false : true;
            }
            
            String manualFlag = getProperty("manual_output");
            if (manualFlag == null) {
                Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.WARNING, "manual_output is not set, will put the file into automatically generated dir");
            } else {
                this.manualOutput = manualFlag.isEmpty() || manualFlag.equalsIgnoreCase("false") ? false : true;
            }

            // Register input files and set up local files
            this.localInputNormalFiles = new String[this.normal.length];
            this.localInputTumorFiles  = new String[this.tumor.length];
            this.normalBases           = new String[this.normal.length];
            this.tumorBases            = new String[this.tumor.length];
            this.bicseqInterval        = Integer.valueOf(getOptionalProperty("biqseq_interval", BICSEQ_I_DEFAULT));
            this.bicseqSpread          = Integer.valueOf(getOptionalProperty("biqseq_spread",   BICSEQ_S_DEFAULT));
                   

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
            Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.SEVERE, null, ex);
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
            Logger.getLogger(BicSeqWorkflow.class.getName()).log(Level.WARNING, null, e);
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
               if (this.templateType.equals(WG)) {
                 // LAUNCH BICseq
                 launchBicSeq(this.localInputNormalFiles[n],
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
     * BICseq configuring/launching
     */
    private void launchBicSeq(String inputNormal, String inputTumor, int id, List<Job> parents) {

        // Job convertJob and create configFile
        Job convertJob = this.getWorkflow().createBashJob("bicseq_prepare");
            
        String configFile = "bicseq_config." + id + ".conf";
        convertJob.setCommand(getWorkflowBaseDir() + "/dependencies/configureBICseq.pl"
                            + " --input-normal " + inputNormal
                            + " --input-tumor "  + inputTumor
                            + " --outdir " + this.dataDir
                            + " --config-file " + configFile
                            + " --samtools " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                            + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/SAMgetUnique/samtools-0.1.7a_getUnique-0.1.1/samtools");
        convertJob.setMaxMemory("4000");
        if (parents != null) {
            for (Job p : parents) {
                convertJob.addParent(p);
            }
        }
        Log.stdout("Created BICseq convert Job");
        
        
        // Launch BICSeq, provision results
        // PERL_pipeline/BICseq_1.1.2/BIC-seq/BIC-seq.pl --I 150,20 /u/pruzanov/Data/CNVtools/BICseq/test1.config /scratch2/users/pruzanov/Data/CNVTOOLS/BIC-seq.hn.test1 \"InitialTest\"
        String resultDir = "BICseq_out_" + id + "/";
        Job launchJob = this.getWorkflow().createBashJob("bicseq_launch");
        String resultID = BICSEQ_PREFIX + this.makeBasename(inputNormal, ".bam") + ".vs." 
                                        + this.makeBasename(inputTumor,  ".bam");
        
        launchJob.setCommand("module load " + this.rModule + ";"
                           + getWorkflowBaseDir() + "/dependencies/launchBICseq.pl"
                           + " --config-file " + this.dataDir + configFile
                           + " --outdir " + this.dataDir + resultDir
                           + " --bicseq-interval " + this.bicseqInterval
                           + " --bicseq-spread "   + this.bicseqSpread
                           + " --result-id " + resultID
                           + " --bicseq " + getWorkflowBaseDir() + "/bin/BICseq-" + this.bicseqVersion
                           + "/PERL_pipeline/BICseq_" + this.bicseqVersion + "/BIC-seq/BIC-seq.pl");
        launchJob.setMaxMemory("6000");
        launchJob.addParent(convertJob);
        Log.stdout("Created BICseq launch Job");
        
        // Provision files normal.vs.tumor.bicseg, normal.vs.tumor.png, normal.vs.tumor.wig
        SqwFile bicseqSegFile = createOutputFile(this.dataDir + resultDir + resultID + ".bicseg", "text/plain", this.manualOutput);
        bicseqSegFile.setSkipIfMissing(skipFlag);
        bicseqSegFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqSegFile);
        
        SqwFile bicseqPngFile = createOutputFile(this.dataDir + resultDir + resultID + ".png",    "image/png",  this.manualOutput);
        bicseqPngFile.setSkipIfMissing(skipFlag);
        bicseqPngFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqPngFile);
        
        SqwFile bicseqWigFile = createOutputFile(this.dataDir + resultDir + resultID + ".wig",    "text/plain", this.manualOutput);
        bicseqWigFile.setSkipIfMissing(skipFlag);
        bicseqWigFile.getAnnotations().put("variation_calling_algorithm", "BICseq " + this.bicseqVersion);
        launchJob.addFile(bicseqWigFile);
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
