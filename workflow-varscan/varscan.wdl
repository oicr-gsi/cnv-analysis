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
call runVarscan { input: inputPileup = makePileups.pileup, sampleID = sampleID }

# Run post-processing job if we have results from runVarscan
Array[File] cNumberFile = select_all([runVarscan.resultFile])
if (length(cNumberFile) == 1) {
    call smoothData{input: copyNumberFile = select_first([runVarscan.resultFile])}
} 


meta {
  author: "Peter Ruzanov"
  email: "peter.ruzanov@oicr.on.ca"
  description: "Varscan 2.0"
}

output {
 File log = runVarscan.logFile
 File? resultFile = smoothData.filteredData
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
 Int? jobMemory = 18
}

parameter_meta {
  inputNormal: "input .bam file for normal tissue"
  inputTumor: "input .bam file for tumor tissue"
  refFasta: "Reference fasta file, path depends on the respective module"
  modules: "required modules"
  jobMemory: "memory for this job, in Gb"
}

command <<<
 ~{samtools} mpileup -q 1 -f ~{refFasta} ~{inputNormal} ~{inputTumor} | awk -F "\t" '$4 > 0 && $7 > 0' > normtumor_sorted.pileup 
>>>

runtime {
 memory: "~{jobMemory} GB"
 modules: "~{modules}"
}

output {
 File pileup = "normtumor_sorted.pileup"
}
}



# ==========================================
#  configure and run Varscan
# ==========================================
task runVarscan {
input {
  File inputPileup
  String? sampleID ="VARSCAN"
  # Parameters
  Float? pValue = 0.05
  Float? delThreshold = 0.5 
  Float? ampThreshold = 0.5 
  Int? minRegionSize = 10
  Int? recenterUp = 0
  Int? recenterDown = 0 
  Int? jobMemory  = 20
  Int? javaMemory = 6
  String? logFile = "VARSCAN.log"
  String? varScan = "$VARSCAN_ROOT/VarScan.jar"
  String? modules = "varscan/2.4.2 java/8"
  Int timeout = 20

}

parameter_meta {
 inputPileup: "Input .pileup file for analysis"
 sampleID: "This is used as a prefix for output files"
 pValue: "p-value for cnv calling, default is 0.05"
 delThreshold: "Upper bound for log ratio to call deletion, default is 0.5"
 ampThreshold: "Lower bound for log ratio to call amplification, default is 0.5"
 minRegionSize: "Minimum size (in bases) for a region to be counted"
 recenterUp: "Recenter data around an adjusted baseline"
 recenterDown: "Recenter data around an adjusted baseline"
 jobMemory: "Memory in Gb for this job"
 javaMemory: "Memory in Gb for Java"
 logFile: "File for logging Varscan messages"
 varScan: "path to varscan .jar file"
 modules: "Names and versions of modules"
 timeout: "Timeout in hours, needed to override imposed limits"
}

command <<<
 module load ~{modules} 2>/dev/null
 python<<CODE
 import os
 import re
 varscan = os.path.expandvars("~{varScan}")
 varscanCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " copynumber ~{inputPileup} ~{sampleID} -mpileup 1 --p-value ~{pValue}"
 cvg = 0
 resultsOk = False
 f = open("VARSCAN.log", "w+")
 f.write('[Varscan log]' + '\n')

 while not resultsOk:
     run_command = varscanCommand
     if cvg != 0:
         run_command += " --min-coverage " + str(cvg)
     else:
         cvg = 19
     cvg -= 4
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
  File logFile   = "VARSCAN.log"
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
 String? min_coverage  = ""
 String? del_threshold = ""
 String? min_region_size = ""
 String? recenter_up = ""
 String? recenter_down = ""
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
 module load ~{modules} 2>/dev/null
 python<<CODE
 import os
 varscan = os.path.expandvars("~{varScan}")
 filterCommand = "java -Xmx~{javaMemory}G -jar " + varscan + " copyCaller ~{copyNumberFile} --output-file ~{sampleID}.copynumber.filtered"

 if "~{min_coverage}":
    filterCommand += " --min-coverage ~{min_coverage}"
 if "~{del_threshold}":
    filterCommand += " --del-threshold ~{del_threshold}"
 if "~{min_region_size}":
    filterCommand += " --min-region-size ~{min_region_size}"
 if "~{recenter_up}":
    filterCommand += " --recenter-up ~{recenter_up}"
 if "~{recenter_down}":
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
