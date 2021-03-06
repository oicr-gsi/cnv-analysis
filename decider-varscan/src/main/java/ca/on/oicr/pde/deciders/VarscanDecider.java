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
public class VarscanDecider extends OicrDecider {
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
    private String forceCrosscheck = "true";
    private String do_sort         = "false";
    private String rmodule         = "R/3.2.1-deb8";
    private String varscanPvalueThreshold;
    private String varscanJavaXmx;
    
    //Additional Parameters for VarScan:
    private String varscanMinCoverage;
    private String varscanDelCoverage;
    private String varscanMinRegion;
    private String varscanRecenterUp;
    private String varscanRecenterDown;

    private final static String BAM_METATYPE = "application/bam";
    private final static String WG           = "WG";
    private final static String EX           = "EX";
    private String tumorType;
    private List<String> duplicates;
    private final static String PVALUE         = "0.05";
    private static final String VARSCAN_JAVA_MEM = "4";
    
    public VarscanDecider() {
        super();
        fileSwaToSmall  = new HashMap<String, BeSmall>();
        parser.acceptsAll(Arrays.asList("ini-file"), "Optional: the location of the INI file.").withRequiredArg();
        parser.accepts("manual-output","Optional*. Set the manual output "
                + "either to true or false").withRequiredArg();
        parser.accepts("template-type","Required. Set the template type to limit the workflow run "
                + "so that it runs on data only of this template type").withRequiredArg();
        parser.accepts("r-module","Optional. Set the R module to load in order to run FREEC scripts ").withRequiredArg();
        parser.accepts("force-crosscheck","Optional. Set the crosscheck to true or false "
                + "when running the workflow, the default is true").withRequiredArg();
        parser.accepts("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./").withRequiredArg();
        parser.accepts("output-folder", "Optional: the name of the folder to put the output into relative to "
	        + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results").withRequiredArg();
        parser.accepts("varscan-min-coverage","Optional. VarScan filtering parameter, see Varscan Manual").withRequiredArg();
        parser.accepts("varscan-del-coverage","Optional. VarScan filtering parameter, see Varscan Manual").withRequiredArg();
        parser.accepts("varscan-min-region","Optional. VarScan filtering parameter, see Varscan Manual").withRequiredArg();
        parser.accepts("varscan-recenter-up","Optional. VarScan filtering parameter, see Varscan Manual").withRequiredArg();
        parser.accepts("varscan-recenter-down","Optional. VarScan filtering parameter, see Varscan Manual").withRequiredArg();
        parser.accepts("queue", "Optional: Set the queue (Default: not set)").withRequiredArg();
        parser.accepts("tumor-type", "Optional: Set tumor tissue type to something other than primary tumor (P), i.e. X . Default: Not set (All)").withRequiredArg();
        parser.accepts("do-sort", "Optional: Set the flag (true or false) to indicate if need to sort bam files. Default: false").withRequiredArg();
        parser.accepts("skip-missing-files","Optional. Set the flag for skipping non-existing files to true or false "
                + "when running the workflow, the default is true").withRequiredArg();
        parser.accepts("varscan-pvalue", "Optional: Set the threshold p-value for Varscan variant calls (0.05 is the default)").withRequiredArg();
        parser.accepts("varscan-java-xmx", "Optional: Set the memory heap in Gigabytes for Varscan java").withRequiredArg();
        parser.accepts("verbose", "Optional: Enable verbose Logging").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        Log.debug("INIT");
	this.setMetaType(Arrays.asList(BAM_METATYPE));
        this.setGroupingStrategy(Header.FILE_SWA);
                
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
                Log.error("--template-type requires an argument, WG or EX");
                rv.setExitStatus(ReturnValue.INVALIDARGUMENT);
                return rv;
            } else {
                this.templateTypeFilter = options.valueOf("template-type").toString();
                if (!this.templateTypeFilter.equals(WG) && !this.templateTypeFilter.equals(EX)) {
                    Log.stderr("NOTE THAT ONLY EX or WG template-type SUPPORTED, WE CANNOT GUARANTEE MEANINGFUL RESULTS WITH OTHER TEMPLATE TYPES");
                }
                this.templateType       = this.templateTypeFilter;
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
        
        if (this.options.has("r-module")) {
            this.rmodule = options.valueOf("r-module").toString();
            Log.debug("Setting R module parameter, default is  R/3.2.1-deb8 and needs to be changed only in special cases");
	}
         
        if (this.options.has("manual-output")) {
            this.manual_output = options.valueOf("manual_output").toString();
            Log.debug("Setting manual output, default is false and needs to be set only in special cases");
	}
        
        
        if (this.options.has("tumor-type")) {
            this.tumorType = options.valueOf("tumor-type").toString();
            Log.debug("Setting tumor type to " + this.tumorType +  " as requested");
	}
        
        if (this.options.has("force-crosscheck")) {
            String crosscheck = options.valueOf("force-crosscheck").toString();
            if (!crosscheck.isEmpty()) {
              this.forceCrosscheck = crosscheck.equalsIgnoreCase("true") ? "true" : "false";
              Log.debug("Setting force crosscheck to " + this.forceCrosscheck);
            }
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
        
        if (options.has("varscan-java-xmx")) {
            this.varscanJavaXmx = options.valueOf("varscan-java-xmx").toString();
        } else {
            this.varscanJavaXmx = VARSCAN_JAVA_MEM;
        }
        
        if (options.has("varscan-pvalue")) {
            this.varscanPvalueThreshold = options.valueOf("varscan-pvalue").toString();
        } else {
            this.varscanPvalueThreshold = PVALUE;
        }
               
        if (options.has("varscan-min-coverage")) {
            this.varscanMinCoverage = options.valueOf("varscan-min-coverage").toString();
        }
        
        if (options.has("varscan-del-coverage")) {
            this.varscanDelCoverage = options.valueOf("varscan-del-coverage").toString();
        }
        
        if (options.has("varscan-min-region")) {
            this.varscanMinRegion = options.valueOf("varscan-min-region").toString();
        }
        
        if (options.has("varscan-recenter-up")) {
            this.varscanRecenterUp = options.valueOf("varscan-recenter-up").toString();
        }
        
        if (options.has("varscan-recenter-down")) {
            this.varscanRecenterDown = options.valueOf("varscan-recenter-down").toString();
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
        int countNorm = 0;
        
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
                    countNorm += 1;
                } else if (!tt.isEmpty()) {
                    haveTumr = true;
                }
            }
        }
        if (countNorm == 1 && haveTumr) {
         return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
        } 
        if (countNorm > 1) {
            Log.error("Multiple Normals detected, WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }
        String absent = haveNorm ? "Tumor" : "Normal";
        Log.error("Data for " + absent + " tissue are not available, WON'T RUN");
        return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
    }

    @Override
    protected boolean checkFileDetails(ReturnValue returnValue, FileMetadata fm) {
        Log.debug("CHECK FILE DETAILS:" + fm);
        String currentTtype    = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle()        + "geo_library_source_template_type");
        String targetResequencingType = returnValue.getAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_targeted_resequencing");
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
        
        //reset test mode
        if (!this.options.has("test")) {
            this.setTest(false);
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
        iniFileMap.put("R_module", this.rmodule);

        iniFileMap.put("manual_output",  this.manual_output);
        iniFileMap.put("force_crosscheck",  this.forceCrosscheck);
        iniFileMap.put("skip_missing_files", this.skipMissing);
        iniFileMap.put("do_sort", this.do_sort);
        if (!this.varscanJavaXmx.isEmpty()) {
            iniFileMap.put("varscan_java_xmx", this.varscanJavaXmx);
        }
        if (!this.varscanPvalueThreshold.isEmpty()) {
            iniFileMap.put("varscan_pvalue", this.varscanPvalueThreshold);
        }
        if (!this.varscanMinCoverage.isEmpty()) {
            iniFileMap.put("varscan_min_coverage", this.varscanMinCoverage);
        }
        if (!this.varscanDelCoverage.isEmpty()) {
            iniFileMap.put("varscan_del_coverage", this.varscanDelCoverage);
        }
        if (!this.varscanMinRegion.isEmpty()) {
            iniFileMap.put("varscan_min_region", this.varscanMinRegion);
        }
        if (!this.varscanRecenterUp.isEmpty()) {
            iniFileMap.put("varscan_recenter_up", this.varscanRecenterUp);
        }
        if (!this.varscanRecenterDown.isEmpty()) {
            iniFileMap.put("varscan_recenter_down", this.varscanRecenterDown);
        }
        
        //Note that we can use group_id, group_description and external_name for tumor bams only
        if (null != groupIds && groupIds.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("group_id", groupIds.toString());
        } else {
          iniFileMap.put("group_id", "NA");    
        }

        if (null != groupDescription && groupDescription.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("group_id_description", groupDescription.toString());
        } else {
          iniFileMap.put("group_id_description", "NA");    
        }
        
        if (null != tubeId && tubeId.length() != 0 && !groupIds.toString().contains("NA")) {
          iniFileMap.put("external_name", tubeId.toString());
        } else {
          iniFileMap.put("external_name", "NA");    
        }
        
        return iniFileMap;
    }
    
   
   public static void main(String args[]){
 
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(VarscanDecider.class.getCanonicalName());
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
