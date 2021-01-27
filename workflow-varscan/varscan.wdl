version 1.0

workflow varscan {
input {
    # Normally we need only tumor bam, normal bam may be used when available
    File inputTumor
    File inputNormal
    File inputTumorIndex
    File inputNormalIndex
    String outputFileNamePrefix = ""
    String bedIntervalsPath = ""
    Array[String] chromRegions = ["chr1:1-249250621","chr2:1-243199373","chr3:1-198022430","chr4:1-191154276","chr5:1-180915260","chr6:1-171115067","chr7:1-159138663","chr8:1-146364022","chr9:1-141213431","chr10:1-135534747","chr11:1-135006516","chr12:1-133851895","chr13:1-115169878","chr14:1-107349540","chr15:1-102531392","chr16:1-90354753","chr17:1-81195210","chr18:1-78077248","chr19:1-59128983","chr20:1-63025520","chr21:1-48129895","chr22:1-51304566","chrX:1-155270560","chrY:1-59373566","chrM:1-16571"]
}

call expandRegions { input: bedPath = bedIntervalsPath }

String sampleID = if outputFileNamePrefix=="" then basename(inputTumor, ".bam") else outputFileNamePrefix
Array[String] splitRegions = if bedIntervalsPath != "" then expandRegions.regions else chromRegions

# Produce pileups
scatter ( r in splitRegions )   {
  call makePileups { input: inputTumor = inputTumor, inputTumorIndex = inputTumorIndex, inputNormal = inputNormal, inputNormalIndex = inputNormalIndex, region = r }
}

# Configure and run Varscan
scatter( p in makePileups.pileup) {
  call runVarscanCNV { input: inputPileup = p, sampleID = sampleID }
  call runVarscanSNV as getSnvNative { input: inputPileup = p, sampleID = sampleID }
  call runVarscanSNV as getSnvVcf { input: inputPileup = p, sampleID = sampleID, outputVcf = 1 }
}

# Merge tasks
call mergeVariantsNative as mergeCNV { input: filePaths = select_all(runVarscanCNV.resultFile), outputFile = sampleID, outputExtension = "copynumber" }
call mergeVariantsNative as mergeSNP { input: filePaths = select_all(getSnvNative.snpFile), outputFile = sampleID, outputExtension = "snp" }
call mergeVariantsNative as mergeIND { input: filePaths = select_all(getSnvNative.indelFile), outputFile = sampleID, outputExtension = "indel" }
call mergeVariantsVcf as mergeSNPvcf { input: filePaths = select_all(getSnvVcf.snpVcfFile), outputSuffix = "snp", outputFile = sampleID }
call mergeVariantsVcf as mergeINDvcf { input: filePaths = select_all(getSnvVcf.indelVcfFile), outputSuffix = "indel", outputFile = sampleID }

# Run post-processing job if we have results from runVarscanCNV
Array[File] cNumberFile = select_all([mergeCNV.mergedVariants])
if (length(cNumberFile) == 1) {
    call smoothData{input: copyNumberFile = select_first([mergeCNV.mergedVariants]), sampleID = sampleID}
}

meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Varscan 2.2, workflow for calling SNVs and CVs"

  dependencies: [
      {
        name: "picard/2.21.2",
        url: "https://broadinstitute.github.io/picard"
      },
      {
        name: "varscan/2.4.2",
        url: "http://varscan.sourceforge.net"
      },
      {
        name: "samtools/0.1.19",
        url: "http://www.htslib.org/"
      },
      {
        name: "rstats/3.6",
        url: "http://cran.utstat.utoronto.ca/src/base/R-3/R-3.6.1.tar.gz"
      }
    ]
    
    output_meta: {
      resultCnvFile: "file with CNV calls, smoothed",
      resultSnpFile: "file with SNPs, native varscan format",
      resultIndelFile: "file with Indel calls, native varscan format",
      resultSnpVcfFile: "file with SNPs, vcf format",
      resultIndelVcfFile: "file with Indels, vcf format"
    }
}

parameter_meta {
  inputTumor: "input .bam file for tumor sample"
  inputNormal: "input .bam file for normal sample"
  inputTumorIndex: "input .bai file for tumor sample"
  inputNormalIndex: "input .bai file for normal sample"
  outputFileNamePrefix: "Output file(s) prefix"
  bedIntervalsPath: "Path to a .bed file used for targeted variant calling"
  chromRegions: "Regions used for scattering tasks, need to be assembly-specific"
}

output {
 File? resultCnvFile      = smoothData.filteredData
 File? resultSnpFile      = mergeSNP.mergedVariants
 File? resultIndelFile    = mergeIND.mergedVariants
 File? resultSnpVcfFile   = mergeSNPvcf.mergedVcf
 File? resultIndelVcfFile = mergeINDvcf.mergedVcf
}

}

# =======================================================
# Read bed file, return a string with regions for mpileup
# =======================================================
task expandRegions {
input {
 String bedPath = ""
 Int jobMemory = 4
}

parameter_meta {
  bedPath: "Optional path to a bed file with intervals"
  jobMemory: "Memory for this task in GB"
}

command <<<
 python <<CODE
 import os
 if os.path.exists("~{bedPath}"):
    with open("~{bedPath}") as f:
        for line in f:
           line = line.rstrip()
           tmp = line.split("\t")
           r = " " + tmp[0] + ":" + tmp[1] + "-" + tmp[2]
           print(r)
    f.close()
 CODE
>>>

runtime {
 memory:  "~{jobMemory} GB"
}

output {
 Array[String] regions = read_lines(stdout()) 
}
}

# ==========================================
#  produce pileup with samtools
# ==========================================
task makePileups {
input {
 File inputNormal
 File inputTumor
 File inputTumorIndex
 File inputNormalIndex
 String refFasta = "$HG19_ROOT/hg19_random.fa"
 String modules  = "samtools/0.1.19 hg19/p13"
 String samtools = "$SAMTOOLS_ROOT/bin/samtools"
 String region 
 Int jobMemory   = 18
 Int timeout     = 40
}

parameter_meta {
  inputNormal: "input .bam file for normal tissue"
  inputNormalIndex: ".bai index file for normal tissue"
  inputTumor: "input .bam file for tumor tissue"
  inputTumorIndex: ".bai index file for tumor tissue"
  refFasta: "Reference fasta file, path depends on the respective module"
  modules: "required modules"
  samtools: "path to samtools"
  region: "Region in a form of chrX:12000-12500 for mpileup command"
  jobMemory: "memory for this job, in Gb"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 set -euxo pipefail
 ~{samtools} mpileup -q 1 -r ~{region} -f ~{refFasta} ~{inputNormal} ~{inputTumor} | awk -F "\t" '$4 > 0 && $7 > 0' | gzip -c > normtumor_sorted.pileup.gz 
>>>

runtime {
 memory: "~{jobMemory} GB"
 modules: "~{modules}"
 timeout: "~{timeout}"
}

output {
 File pileup = "normtumor_sorted.pileup.gz"
}
}

#=============================================================
# Task for concatenating CNV and SNV variants
#=============================================================
task mergeVariantsNative {
input {
 Array[File] filePaths
 String outputFile = "concatenated_variants"
 String outputExtension = "csv"
 Int jobMemory = 6
 Int timeout   = 10
}

parameter_meta {
  filePaths: "Array of pileup files to concatenate"
  jobMemory: "memory in GB for this job"
  outputExtension: "Extension of the output file"
  outputFile: "Name of the output file"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 set -euo pipefail
 head -n 1 ~{filePaths[0]} > "~{outputFile}.~{outputExtension}"
 cat ~{sep=' ' filePaths} | sort -V -k 1,2 | grep -v ^chrom | grep -v ^chrM >> "~{outputFile}.~{outputExtension}"
 cat ~{sep=' ' filePaths} | awk '{if($1 == "chrM"){print $0}}' | sort -V -k 1,2 >> "~{outputFile}.~{outputExtension}"
 if [ ! -s ~{outputFile}.~{outputExtension} ] ; then
  rm ~{outputFile}.~{outputExtension}
 fi
>>>

runtime {
 memory: "~{jobMemory} GB"
 timeout: "~{timeout}"
}

output {
  File? mergedVariants = "~{outputFile}.~{outputExtension}"
}
}


#=============================================================
# Task for concatenating CNV and SNV variants
#=============================================================
task mergeVariantsVcf {
input {
 Array[File] filePaths
 String outputFile = "concatenated_vcf"
 String outputSuffix = "snp"
 String modules = "picard/2.21.2 hg19/p13"
 String seqDictionary = "$HG19_ROOT/hg19_random.dict"
 Int jobMemory = 12
 Int javaMemory = 8
 Int timeout   = 10
}

parameter_meta {
  filePaths: "Array of pileup files to concatenate"
  jobMemory: "memory in GB for this job"
  javaMemory: "memory in GB for java VM"
  outputFile: "Name of the output file"
  outputSuffix: "Suffix to use for an output file: snp or indel"
  seqDictionary: ".dict file for the reference in use"
  modules: "modules needed for this task"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command<<<
 set -euxo pipefail
 unset _JAVA_OPTIONS
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar SortVcf I=~{sep=' I=' filePaths} SD=~{seqDictionary} O=~{outputFile}.~{outputSuffix}.vcf
>>>

runtime {
 modules: "~{modules}"
 memory: "~{jobMemory} GB"
 timeout: "~{timeout}"
}

output {
  File? mergedVcf = "~{outputFile}.~{outputSuffix}.vcf"
}
}

# ==========================================
#  configure and run Varscan in SNV mode
# ==========================================
task runVarscanSNV {
input {
  File inputPileup
  String sampleID ="VARSCAN"
  Float pValue = 0.05
  Int jobMemory  = 20
  Int javaMemory = 6
  Int minCoverage = 8
  Int minCoverageNormal = 8
  Int minCoverageTumor = 6
  Float minVarFreq = 0.1
  Float minFreqForHom = 0.75
  Float normalPurity = 1.0
  Float tumorPurity = 1.0
  Float pValueHet = 0.99
  Int strandFilter = 0
  Int validation = 0
  Int outputVcf = 0
  String logFile = "VARSCAN_SNV.log"
  String varScan = "$VARSCAN_ROOT/VarScan.jar"
  String modules = "varscan/2.4.2 java/8"
  Int timeout = 40
}

parameter_meta {
 inputPileup: "Input .pileup file for analysis"
 sampleID: "This is used as a prefix for output files"
 pValue: "somatic p-value for SNV calling, default is 0.05"
 minCoverage: "Minimum coverage in normal and tumor to call variant [8]"
 minCoverageNormal: "Minimum coverage in normal to call somatic [8]"
 minCoverageTumor: "Minimum coverage in tumor to call somatic [6]"
 minVarFreq: "Minimum variant frequency to call a heterozygote [0.10]"
 minFreqForHom: "Minimum frequency to call homozygote [0.75]"
 normalPurity: "Estimated purity (non-tumor content) of normal sample [1.00]"
 tumorPurity: "Estimated purity (tumor content) of normal sample [1.00]"
 pValueHet: "p-value threshold to call a heterozygote [0.99]"
 strandFilter: "If set to 1, removes variants with >90% strand bias"
 validation: "If set to 1, outputs all compared positions even if non-variant"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 outputVcf: "Flag that when set to 1 indicates that we need results in vcf format"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 unset _JAVA_OPTIONS
 set -euxo pipefail
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "zcat ~{inputPileup} | java -Xmx~{javaMemory}G -jar " + varscan + " somatic -mpileup 1 --somatic-p-value ~{pValue}"

 if "~{minCoverageNormal}" != "8":
    varscanCommand += " --min-coverage-normal ~{minCoverageNormal}"
 if "~{minCoverageTumor}" != "6":
    varscanCommand += " --min-coverage-tumor ~{minCoverageTumor}"
 if "~{minVarFreq}" != "0.1":
    varscanCommand += " --min-var-freq ~{minVarFreq}"
 if "~{minFreqForHom}" != "0.75":
    varscanCommand += " --min-freq-for-hom ~{minFreqForHom}"
 if "~{normalPurity}" != "1.0":
    varscanCommand += " --normal-purity ~{normalPurity}"
 if "~{tumorPurity}" != "1.0":
    varscanCommand += " --tumor-purity ~{tumorPurity}"
 if "~{pValueHet}" != "0.99":
    varscanCommand += " --p-value ~{pValueHet}"
 if "~{strandFilter}" != "0":
    varscanCommand += " --strand-filter 1"
 if "~{validation}" != "0":
    varscanCommand += " --validation 1"
 if "~{outputVcf}" != "0":
    varscanCommand += " --output-vcf 1"
    varscanCommand += " --output-snp ~{sampleID}.snp.vcf --output-indel ~{sampleID}.indel.vcf"
 else:
    varscanCommand += " --output-snp ~{sampleID}.snp --output-indel ~{sampleID}.indel"

 cvg = ~{minCoverage}
 resultsOk = False
 f = open("~{logFile}", "w+")
 f.write('[Varscan log]' + '\n')

 while not resultsOk:
     run_command = varscanCommand
     run_command += " --min-coverage " + str(cvg)
     cvg -= 2
     f.write('[' + run_command + ']\n')
     res_string = os.popen(run_command + " 2>&1").read()
     f.write(res_string)
     m = re.search(r'(\d+)', str(res_string))
     m = re.search(r'(\d+) had sufficient coverage', res_string)
     if not m or (m and m.group(1) == '0' and cvg <= 2):
         f.write('Unable to run Varscan even with min-coverage set to ' + str(cvg) + ' aborting...\n')
         break
     elif m and m.group(1) != '0':
         resultsOk = True
         f.write('Success!\n')
     else:
         f.write('Coverage threshold too high, trying min-coverage ' + str(cvg) + '...\n')

 if not resultsOk:
     f.write('Varscan failed\n')
 f.close()
 CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File? snpFile   = "~{sampleID}.snp"
  File? indelFile = "~{sampleID}.indel"
  File? snpVcfFile = "~{sampleID}.snp.vcf"
  File? indelVcfFile = "~{sampleID}.indel.vcf"
}
}

# ==========================================
#  configure and run Varscan in CNV mode
# ==========================================
task runVarscanCNV {
input {
  File inputPileup
  String sampleID ="VARSCAN"
  Float pValue = 0.05
  Int jobMemory  = 20
  Int javaMemory = 6
  String logFile = "VARSCAN_CNV.log"
  String varScan = "$VARSCAN_ROOT/VarScan.jar"
  String modules = "varscan/2.4.2 java/8"
  Int timeout = 40
}

parameter_meta {
 inputPileup: "Input .pileup file for analysis"
 sampleID: "This is used as a prefix for output files"
 pValue: "p-value for cnv calling, default is 0.05"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 unset _JAVA_OPTIONS
 set -euxo pipefail
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "zcat ~{inputPileup} | java -Xmx~{javaMemory}G -jar " + varscan + " copynumber --output-file ~{sampleID} -mpileup 1 --p-value ~{pValue}"
 cvg = 0
 resultsOk = False
 f = open("~{logFile}", "w+")
 f.write('[Varscan log]' + '\n')

 while not resultsOk:
     run_command = varscanCommand
     if cvg == 0:
         cvg = 19
     cvg -= 4
     run_command += " --min-coverage " + str(cvg)
     f.write('[' + run_command + ']\n')
     res_string = os.popen(run_command + " 2>&1").read()
     m = re.search(r'(\d+)', str(res_string))
     m = re.search(r'(\d+) had sufficient coverage', res_string)
     if not m or (m and m.group(1) == '0' and cvg <= 2):
         f.write('Unable to run Varscan even with min-coverage set to ' + str(cvg) + ' aborting...\n')
         break
     elif m and m.group(1) != '0':
         resultsOk = True
         f.write('Success!\n')
     else:
         f.write('Coverage threshold too high, trying min-coverage ' + str(cvg) + '...\n')

 if not resultsOk:
     f.write('Varscan failed\n')
 f.close()
 CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File? resultFile = "~{sampleID}.copynumber"
}
}

# ====================================================
#  Additional smoothing for Varscan data. Visualization
#  will be implemented after updates to rstats module
#       needed to enable png/jpg rendering
# ======================================================

task smoothData {
input {
 File copyNumberFile
 String varScan = "$VARSCAN_ROOT/VarScan.jar"
 String modules = "varscan/2.4.2 java/8 rstats/3.6"
 Int min_coverage  = 20
 Int max_homdel_coverage = 5
 Int min_tumor_coverage = 10
 Float del_threshold = 0.25
 Float amp_threshold = 0.25
 Int min_region_size = 10
 Int recenter_up = 0
 Int recenter_down = 0
 String sampleID ="VARSCAN"
 Int jobMemory  = 16
 Int javaMemory = 6
}

parameter_meta {
 copyNumberFile: "Output from Varscan"
 varScan: "Path to VarScan jar file"
 modules: "Modules for this job"
 min_coverage: "Fine-tuning parameter for VarScan"
 max_homdel_coverage: "Max coverage form homozygous deletion, default is 5"
 min_tumor_coverage: "Min coverage in tumor sample, default is 10"
 del_threshold: "Fine-tuning parameter for VarScan"
 amp_threshold: "Amplification threshold to report, default is 0.25"
 min_region_size: "Fine-tuning parameter for VarScan"
 recenter_up: "Fine-tuning parameter for VarScan"
 recenter_down: "Fine-tuning parameter for VarScan"
 sampleID: "sample id (used as prefix for result files)"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
}

command <<<
 python<<CODE
 import os
 varscan = os.path.expandvars("~{varScan}")
 filterCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " copyCaller ~{copyNumberFile} --output-file ~{sampleID}.copynumber.filtered"

 if "~{min_coverage}" != "20":
    filterCommand += " --min-coverage ~{min_coverage}"
 if "~{min_tumor_coverage}" != "10":
    filterCommand += " --min-tumor-coverage ~{min_tumor_coverage}"
 if "~{max_homdel_coverage}" != "5":
    filterCommand += "--max-homdel-coverage ~{max_homdel_coverage}"
 if "~{del_threshold}" != "0.25":
    filterCommand += " --del-threshold ~{del_threshold}"
 if "~{amp_threshold}" != "0.25":
    filterCommand += " --amp-threshold ~{amp_threshold}"
 if "~{min_region_size}" != "10":
    filterCommand += " --min-region-size ~{min_region_size}"
 if "~{recenter_up}" != "0":
    filterCommand += " --recenter-up ~{recenter_up}"
 if "~{recenter_down}" != "0":
    filterCommand += " --recenter-down ~{recenter_down}"

 message = os.popen(filterCommand + " 2>&1").read()
 CODE
>>>

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
}

output {
 File? filteredData = "~{sampleID}.copynumber.filtered"
}

}

