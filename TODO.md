# Problems
1. Bwa index fails to add the -p prefix in conda setting.
2. When skip-trimming option is picked the workflow fails due to bad naming practice in the workflow.
3. Conda option did not work before and fixed now, but can be made leaner.
4. Unnecessary R packages exist in the workflow such as DESeq2 and edgeR. Pick one and remove the other.

# TODO
- [x] Fix the bwa index issue.
Issue solved by moving the -p first and fasta file as the last argument.
- [] Fix the naming convention in trim-galore output and bwa input. Try pulling all the names from the metadata file by defining objects.

def create_fastq_channel(LinkedHashMap row) also needs to be adjusted for below to work.
'''
# this is already used to pull the prefix from metadata file
def name = task.ext.prefix ?: "${meta.sample_id}"
# Include fastq file names such as
def file1 = "${meta.file1}"
def file2 = "${meta.file2}"
'''

- [ ] Remove unnecessary packages from the workflow to make conda leaner.
- [ ] Refactor the R code to make it more readable.
- [ ] Consider replacing the Python code with R to reduce workflow dependencies.