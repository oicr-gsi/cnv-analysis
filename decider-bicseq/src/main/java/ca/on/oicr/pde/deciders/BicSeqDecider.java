/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package ca.on.oicr.pde.deciders;

import java.io.File;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.FileMetadata;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 *
 * @author pruzanov@oicr.on.ca
 * 
 * Before running test (i.e 
 *  mvn failsafe:integration-test -DskipITs=false 
 *                                -DwebserviceUrl=http://hsqwstage-www1.hpc.oicr.on.ca:8080/seqware-webservice)
 * issue this command:
 * export _JAVA_OPTIONS="-Xmx3000M"
 */
public class BicSeqDecider extends OicrDecider {
    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Map<String, BeSmall> fileSwaToSmall;

    
    //CNV specific stuff
    private String templateTypeFilter = "";
    private String templateType       = "";
    private String output_prefix      = "./";
    private String queue              = " ";
    private String output_dir      = "seqware-results";
    private String skipMissing     = "true";
    private String manual_output   = "false";
    private String do_sort         = "false";
    private String rmodule         = "R/3.2.1-deb8";
    private String biqseqInterval  = "";
    private String biqseqSpread    = "";
    private static final String BICSEQ_I_DEFAULT         = "150";
    private static final String BICSEQ_S_DEFAULT         = "20";

    private final static String BAM_METATYPE = "application/bam";
    private final static String WG           = "WG";
    private String tumorType;
    private List<String> duplicates;
    
    public BicSeqDecider() {
        super();
        fileSwaToSmall  = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("manual-output","Optional*. Set the manual output "
                + "either to true or false").withRequiredArg();
        parser.accepts("template-type","Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("r-module","Optional. Set the R module to load in order to run HMMcopy scripts ").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
	        + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("biqseq-interval", "Optional: Interval parameter used by BicSeq (Default: 150)").withRequiredArg();
        parser.accepts("biqseq-spread", "Optional: Spread parameter used by BicSeq (Default: 20)").withRequiredArg();       
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("tumor-type", "Optional: Set tumor tissue type to something other than primary tumor (P), i.e. X . Default: Not set (All)").withRequiredArg();
        parser.accepts("do-sort", "Optional: Set the flag (true or false) to indicate if need to sort bam files. Default: false").withRequiredArg();
        parser.accepts("skip-missing-files","Optional. Set the flag for skipping non-existing files to true or false "
                + "when running the workflow, the default is true").withRequiredArg();
        parser.accepts("verbose", "Optional: Enable verbose Logging").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
	this.setMetaType(Arrays.asList(BAM_METATYPE));
                
        ReturnValue rv = super.init();
        rv.setExitStatus(ReturnValue.SUCCESS);
        
	//Group by sample if no other grouping selected
        if (this.options.has("group-by")) {
            Log.error("group-by parameter passed, but this decider does not allow overriding the default grouping (by Donor + Library Type)");
        }
               
        if (this.options.has("queue")) {
            this.queue   = options.valueOf("queue").toString();
	} else {
            this.queue   = " ";
        }
         
        this.templateTypeFilter = WG;
        if (this.options.has("template-type")) {
            if (!options.hasArgument("template-type")) {
                Log.error("--template-type requires an argument");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            } else {
                this.templateTypeFilter = options.valueOf("template-type").toString();
                if (!this.templateTypeFilter.equals(WG)) {
                    Log.stderr("NOTE THAT ONLY WG template-type SUPPORTED");
                    rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                    return rv;
                }
                this.templateType = this.templateTypeFilter;
            }
	}
        
         if (this.options.has("skip-missing-files")) {
            if (options.hasArgument("skip-missing-files")) {
                this.skipMissing = options.valueOf("skip-missing-files").toString();
                if (!this.skipMissing.equals("false")) {
                    this.skipMissing = "true"; // Default is true, so we care only when it is set to false
                }
            } 
	}
        
        if (this.options.has("manual-output")) {
            this.manual_output = options.valueOf("manual_output").toString();
            Log.debug("Setting manual output, default is false and needs to be set only in special cases");
	}
        
        if (this.options.has("r-module")) {
            this.rmodule = options.valueOf("r-module").toString();
            Log.debug("Setting R module parameter, default is  R/3.2.1-deb8 and needs to be changed only in special cases");
	}
        
        if (this.options.has("tumor-type")) {
            this.tumorType = options.valueOf("tumor-type").toString();
            Log.debug("Setting tumor type to " + this.tumorType +  " as requested");
	}
        
        if (this.options.has("biqseq-interval")) {
            this.biqseqInterval = options.valueOf("biqseq-interval").toString();
            Log.debug("Setting BiqSeq interval, default is " + BICSEQ_I_DEFAULT);
	}
        
        if (this.options.has("biqseq-spread")) {
            this.biqseqSpread = options.valueOf("biqseq-spread").toString();
            Log.debug("Setting BiqSeq spread, default is " + BICSEQ_S_DEFAULT);
	}  
       
        if (this.options.has("verbose")) {
            Log.setVerbose(true);
	} 
        
        if (this.options.has("do-sort")) {
            String tempSort = options.valueOf("do-sort").toString();
            if (tempSort.equalsIgnoreCase("false") || tempSort.equalsIgnoreCase("true"))
                this.do_sort = tempSort.toLowerCase();
        }

        if (this.options.has("output-path")) {
             this.output_prefix = options.valueOf("output-path").toString();
              if (!this.output_prefix.endsWith("/")) {
                 this.output_prefix += "/";
              }
        }
        
        if (this.options.has("output-folder")) {
            this.output_dir = options.valueOf("output-folder").toString();
	}
        
        // Warn about using force-run-all (may not be relevant after 1.0.17 release)
        if (options.has("force-run-all")) {
            Log.stderr("Using --force-run-all WILL BREAK THE LOGIC OF THIS DECIDER, USE AT YOUR OWN RISK");
        }

        return rv;
    }

    /**
     * Final check
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     * @return    */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {
        String[] filePaths = commaSeparatedFilePaths.split(",");
        boolean haveNorm = false;
        boolean haveTumr = false;
        
        // Check for duplicate file names and exclude them from analysis
        this.duplicates = detectDuplicates(commaSeparatedFilePaths);
        
        for (String p : filePaths) {
            if (null != this.duplicates && this.duplicates.contains(p)) {
                Log.stderr("File [" + p + "] has a name that cannot be disambiguated in current set, will skip it");
                continue;
            }
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p))
                    continue;
                String tt = bs.getTissueType();
                

                if (!tt.isEmpty() && tt.equals("R")) {
                    haveNorm = true;
                } else if (!tt.isEmpty()) {
                    haveTumr = true;
                }
            }
        }
        if (haveNorm && haveTumr) {
         return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        } 
            
        String absent = haveNorm ? "Tumor" : "Normal";
        Log.error("Data for " + absent + " tissue are not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype    = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle()        + "geo_library_source_template_type");
        String currentTissueType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle()      + "geo_tissue_type" );

        if (null == currentTissueType )
            return false; // we need only those which have their tissue type set
        // Filter the data of a different template type if filter is specified
        if (!this.templateTypeFilter.equalsIgnoreCase(currentTtype))
            return false;
        // Do not process tumor tissues of type that doesn't match set parameter
        if (null != this.tumorType) {
          if (!currentTissueType.equals("R") && !currentTissueType.equals(this.tumorType))
            return false;
        }
        
        if (this.templateType.isEmpty() || !this.templateType.equals(currentTtype)) {
            this.templateType = currentTtype;
        }
                 
        for (FileMetadata fmeta : returnValue.getFiles()) {
            if (!fmeta.getMetaType().equals(BAM_METATYPE))
                continue;
            if (!fmeta.getFilePath().contains("sorted"))
                this.do_sort = "true"; // Force sorting of all files even if only one is unsorted
        }
       
        return super.checkFileDetails(returnValue, fm);
    }

    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        // get files from study
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();
        // Override the supplied group-by value
        for (ReturnValue currentRV : vals) {
            boolean metatypeOK = false;
            
            for (int f = 0; f < currentRV.getFiles().size(); f++) {
               try {
                 if (currentRV.getFiles().get(f).getMetaType().equals(BAM_METATYPE))
                     metatypeOK = true;
               } catch (Exception e) {
                 Log.stderr("Error checking a file");
               }
            }
            
            if (!metatypeOK)
                continue; // Go to the next value
            
            
            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);
            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();
            
            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                iusDeetsToRV.put(fileDeets, currentRV);
            } 
            //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);               
                BeSmall oldSmall  = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate      = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    iusDeetsToRV.put(fileDeets, currentRV);
                } 
            }
        }

        
        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());       
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();

        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();
            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }
        
        return map;
    }

    @Override
    protected String handleGroupByAttribute(String attribute) {
        String a = super.handleGroupByAttribute(attribute);
        BeSmall small = fileSwaToSmall.get(a);
        if (small != null) {
            return small.getGroupByAttribute();
        }
        return attribute;
    }

    @Override
    protected Map<String, String> modifyIniFile(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {

        StringBuilder inputNormFiles  = new StringBuilder();
        StringBuilder inputTumrFiles  = new StringBuilder();
        StringBuilder groupIds        = new StringBuilder();
        String [] filePaths = commaSeparatedFilePaths.split(",");
        StringBuilder tubeId          = new StringBuilder();
        StringBuilder groupDescription= new StringBuilder();
        
        for (String p :  filePaths) {
            if (null != this.duplicates && this.duplicates.contains(p)) {
                Log.stderr("Will not include file [" + p + "] since there is an ambiguity in names that cannot be resolved");
                continue;
            }
            
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p))
                    continue;

                String tt = bs.getTissueType();
                if (!tt.isEmpty() && tt.equals("R")) {
                    if (inputNormFiles.length() != 0) {
                     inputNormFiles.append(",");
                    }
                 inputNormFiles.append(p);
                } else if (!tt.isEmpty()) {
                    if (inputTumrFiles.length() != 0) {
                     inputTumrFiles.append(",");
                     // group_ids recoreded using info from tumor entries, normal files do not have group_ids
                     groupIds.append(",");
                     groupDescription.append(",");
                     tubeId.append(",");
                    }
                 inputTumrFiles.append(p);
                 groupIds.append(bs.getGroupID());
                 groupDescription.append(bs.getGroupDescription());
                 tubeId.append(bs.getTubeId());
                }
            }
        }
        
      
        // Just in case
        // This should handle possible problems with --force-run-all
        if (inputNormFiles.length() == 0 || inputTumrFiles.length() == 0) {
         Log.error("THE DONOR does not have data to run the workflow");
         this.setTest(true);
        }       
        
        Map<String, String> iniFileMap = new TreeMap<String, String>();
        
        iniFileMap.put("input_files_normal", inputNormFiles.toString());
        iniFileMap.put("input_files_tumor",  inputTumrFiles.toString());
        iniFileMap.put("data_dir", "data");
        iniFileMap.put("template_type", this.templateType);
 
	iniFileMap.put("output_prefix",this.output_prefix);
	iniFileMap.put("output_dir", this.output_dir);
        if (!this.queue.isEmpty()) {
         iniFileMap.put("queue", this.queue);
        } else {
         iniFileMap.put("queue", " ");
        }

        iniFileMap.put("manual_output",  this.manual_output);
        iniFileMap.put("skip_missing_files", this.skipMissing);
        iniFileMap.put("do_sort", this.do_sort);
        iniFileMap.put("R_module", this.rmodule);
        if (!this.biqseqInterval.isEmpty())
            iniFileMap.put("biqseq_interval", this.biqseqInterval);
        if (!this.biqseqSpread.isEmpty())
            iniFileMap.put("biqseq_spread", this.biqseqSpread);
        
        
        //Note that we can use group_id, group_description and external_name for tumor bams only
        if (groupIds.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("group_id", groupIds.toString());
        } else {
          iniFileMap.put("group_id", "NA");    
        }

        if (groupDescription.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("group_id_description", groupDescription.toString());
        } else {
          iniFileMap.put("group_id_description", "NA");    
        }
        
        if (tubeId.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("external_name", tubeId.toString());
        } else {
          iniFileMap.put("external_name", "NA");    
        }
        
        return iniFileMap;
    }
    
   
   public static void main(String args[]){
 
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(BicSeqDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));
         
    }
   
   private class BeSmall {

        private Date   date = null;
        private String iusDetails = null;
        private String groupByAttribute = null;
        private String tissueType = null;
        private String path = null;
        private String tubeID = null;
        private String groupID = null;
        private String groupDescription = null;
        
        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            tissueType = fa.getLimsValue(Lims.TISSUE_TYPE);
            tubeID     = fa.getLimsValue(Lims.TUBE_ID);
            if (null == tubeID || tubeID.isEmpty()) {
                tubeID = "NA";
            }
            groupID    = fa.getLimsValue(Lims.GROUP_ID);
            if (null == groupID || groupID.isEmpty()) {
                groupID = "NA";
            }
            groupDescription = fa.getLimsValue(Lims.GROUP_DESC);
            if (null == groupDescription || groupDescription.isEmpty()) {
                groupDescription = "NA";
            }
            groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);           
            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }
        
        public String getTissueType() {
            return tissueType;
        }
        
        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }
        
        public String getTubeId () {
            return tubeID;
        }
        
        public String getGroupID() {
            return groupID;
        }

        public String getGroupDescription() {
            return groupDescription;
        }

        public void setPath(String path) {
            this.path = path;
        }
    }
   public static boolean fileExistsAndIsAccessible(String filePath) {

        File file = new File(filePath);
        return (file.exists() && file.canRead() && file.isFile());

    }
   
   public static List<String> detectDuplicates(String commaSeparatedFilePaths) {
       
       String [] filePaths = commaSeparatedFilePaths.split(",");
       List<String> list    = new ArrayList<String>();
       List<String> checker = new ArrayList<String>();
       
       for (String path : filePaths) {
           String baseName = makeBasename(path, ".bam");
           
           if (checker.contains(baseName) && !list.contains(path)) {
               list.add(path);
           } else {
               checker.add(baseName);
           }
       }
       
       return list.isEmpty() ? null : list;
       
   }
   
    /**
     * Utility function
     * 
     * @param path
     * @param extension
     * @return 
     */
    public static String makeBasename(String path, String extension) {
        return path.substring(path.lastIndexOf("/") + 1, path.lastIndexOf(extension));
    }
}
