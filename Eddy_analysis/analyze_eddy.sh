
#!/bin/bash

folder=${1}
../fix_files_symbols.sh ${folder}
outfolder="output_$(basename ${folder})"
mkdir -p ${outfolder}
root -l -q plot_aligned_kicks.C\(\"${folder}\"\,\"${outfolder}\"\)
root -l plot_output.C\(\"${outfolder}/output.root\"\)

