#########################################
###  					 Inputs
#########################################
input = ../input/EKNR_Simple_2Cty_Shocks.mat

##########################################
###   					Outputs
##########################################
folders_and_symlinks = ../output/done_create_folders_and_symlinks.txt

pdf = ../output/EKNR_2Cty.pdf

output_all += $(pdf)

#########################################
###   				Work Session
#########################################
all : $(folders_and_symlinks) $(output_all)

$(folders_and_symlinks) : create_folders_and_symlinks.sh
	-bash create_folders_and_symlinks.sh
	echo ""> ../output/done_create_folders_and_symlinks.txt

$(output_all) : $(input) Simple_2Cty.jl Simple_2Cty_Functions.jl
	-julia -q Simple_2Cty.jl
