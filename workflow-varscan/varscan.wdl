version 1.0

workflow varscan {
input {
    # Normally we need only tumor bam, normal bam may be used when available
    File inputTumor
    File inputNormal
    String? outputFileNamePrefix = ""
}

String? sampleID = if outputFileNamePrefix=="" then basename(inputTumor, ".bam") else outputFileNamePrefix

# Produce pileups
call makePileups { input: inputTumor = inputTumor, inputNormal = inputNormal }

# Configure and run Varscan
call runVarscanCNV { input: inputPileup = makePileups.pileup, sampleID = sampleID }

call runVarscanSNV as getSnvNative { input: inputPileup = makePileups.pileup, sampleID = sampleID }
call runVarscanSNV as getSnvVcf { input: inputPileup = makePileups.pileup, sampleID = sampleID, outputVcf = 1 }

# Run post-processing job if we have results from runVarscanCNV
Array[File] cNumberFile = select_all([runVarscanCNV.resultFile])
if (length(cNumberFile) == 1) {
    call smoothData{input: copyNumberFile = select_first([runVarscanCNV.resultFile]), sampleID = sampleID}
}


meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Varscan 2.0"
}

output {
 File? resultCnvFile      = smoothData.filteredData
 File? resultSnpFile      = getSnvNative.snpFile
 File? resultIndelFile    = getSnvNative.indelFile
 File? resultSnpVcfFile   = getSnvVcf.snpVcfFile
 File? resultIndelVcfFile = getSnvVcf.indelVcfFile
}

}

# ==========================================
#  produce pileup with samtools
# ==========================================
task makePileups {
input {
 File inputNormal
 File inputTumor
 String? refFasta = "$HG19_ROOT/hg19_random.fa"
 String? modules  = "samtools/0.1.19 hg19/p13"
 String? samtools = "$SAMTOOLS_ROOT/bin/samtools"
 Int? jobMemory   = 18
 Int timeout      = 20
}

parameter_meta {
  inputNormal: "input .bam file for normal tissue"
  inputTumor: "input .bam file for tumor tissue"
  refFasta: "Reference fasta file, path depends on the respective module"
  modules: "required modules"
  jobMemory: "memory for this job, in Gb"
  timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 ~{samtools} mpileup -q 1 -f ~{refFasta} ~{inputNormal} ~{inputTumor} | awk -F "\t" '$4 > 0 && $7 > 0' > normtumor_sorted.pileup 
>>>

runtime {
 memory: "~{jobMemory} GB"
 modules: "~{modules}"
 timeout: "~{timeout}"
}

output {
 File pileup = "normtumor_sorted.pileup"
}
}

# ==========================================
#  configure and run Varscan in SNV mode
# ==========================================
task runVarscanSNV {
input {
  File inputPileup
  String? sampleID ="VARSCAN"
  Float? pValue = 0.05
  Int? jobMemory  = 20
  Int? javaMemory = 6
  Int? minCoverage       = 8
  Int? minCoverageNormal = 8
  Int? minCoverageTumor  = 6
  Float? minVarFreq        = 0.1
  Float? minFreqForHom     = 0.75
  Float? normalPurity      = 1.0
  Float? tumorPurity       = 1.0
  Float? pValueHet         = 0.99
  Int strandFilter         = 0
  Int validation           = 0
  Int? outputVcf           = 0
  String? logFile          = "VARSCAN_SNV.log"
  String? varScan          = "$VARSCAN_ROOT/VarScan.jar"
  String? modules          = "varscan/2.4.2 java/8"
  Int timeout              = 20
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
 pValueHet: "p-value threshold to call a heterozygote [0.99]"
 strandFilter: "If set to 1, removes variants with >90% strand bias"
 validation: "If set to 1, outputs all compared positions even if non-variant"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 unset _JAVA_OPTIONS
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " somatic ~{inputPileup} ~{sampleID} -mpileup 1 --somatic-p-value ~{pValue}"

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
  String? sampleID ="VARSCAN"
  Float? pValue = 0.05
  Int? jobMemory  = 20
  Int? javaMemory = 6
  String? logFile = "VARSCAN_CNV.log"
  String? varScan = "$VARSCAN_ROOT/VarScan.jar"
  String? modules = "varscan/2.4.2 java/8"
  Int timeout = 20
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
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " copynumber ~{inputPileup} ~{sampleID} -mpileup 1 --p-value ~{pValue}"
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
 String? varScan = "$VARSCAN_ROOT/VarScan.jar"
 String? modules = "varscan/2.4.2 java/8 rstats/3.6"
 String? min_coverage  = "20"
 String? max_homdel_coverage = "5"
 String? min_tumor_coverage = "10"
 String? del_threshold = "0.25"
 String? amp_threshold = "0.25"
 String? min_region_size = "10"
 String? recenter_up = "0"
 String? recenter_down = "0"
 String? sampleID ="VARSCAN"
 Int? jobMemory  = 16
 Int? javaMemory = 6
}

parameter_meta {
 copyNumberFile: "Output from Varscan"
 varScan: "Path to VarScan jar file"
 modules: "Modules for this job"
 min_coverage: "Fine-tuning parameter for VarScan"
 del_threshold: "Fine-tuning parameter for VarScan"
 min_region_size: "Fine-tuning parameter for VarScan"
 recenter_up: "Fine-tuning parameter for VarScan"
 recenter_down: "Fine-tuning parameter for VarScan"
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
