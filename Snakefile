from pathlib import Path
SUBJECTS = sorted(p.name for p in Path("mri_dataset").glob("sub-*"))
SESSIONS = {
  subject: [
    p.name for p in Path(f"mri_dataset/{subject}").glob("ses-*")
    if (p / f"mixed/{subject}_{p.name}_acq-mixed_SE-modulus.nii.gz").exists()
  ]
  for subject in SUBJECTS
}
rule all:
  input:
    [
      f"mri_processed_data/{subject}/concentrations/{subject}_{ses}_acq-mixed_concentration.nii.gz"
      for subject in SUBJECTS
      for ses in SESSIONS[subject]
    ]
  run:
    print("Running pipeline with data from the following subjects and sessions:")
    for subject, sessions in SESSIONS.items():
      print(f"{subject}: {sessions}")


# Snakemake creates output directories, before running the shell command.
# Since FreeSurfer checks for existing subject directories before running
# recon_all with the usual command, I here do and alternative approach,
# where we copy the T1 and FLAIR files into the target directories, and 
# start from there.
rule recon_all:
  input:
    t1w="mri_dataset/{subject}/ses-01/anat/sub-01_ses-01_T1w.nii.gz",
    flair="mri_dataset/{subject}/ses-01/anat/sub-01_ses-01_FLAIR.nii.gz",
  output: # protected to avoid accidental reruns.
    seg=protected("mri_processed_data/freesurfer/{subject}/mri/aparc+aseg.mgz"),
    t1w=protected("mri_processed_data/freesurfer/{subject}/mri/orig/001.mgz"),
    flair=protected("mri_processed_data/freesurfer/{subject}/mri/orig/FLAIRraw.mgz"),
  container:
    "docker://freesurfer/freesurfer:7.4.1"
  shell:
    "mri_convert {input.t1w} {output.t1w} && "  
    "mri_convert --no_scale 1 {input.flair} {output.flair} && "
    "recon-all"
    " -all"
    " -sd mri_processed_data/freesurfer"
    " -s {wildcards.subject}"
    " -parallel"
    " -FLAIRpial"
    " -FLAIR {output.flair}"

rule estimate_T1maps:
  input:
    se="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_SE-modulus.nii.gz",
    ir="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_IR-corrected-real.nii.gz",
    meta="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_meta.json",
  output:
    raw="mri_dataset/derivatives/{subject}/{ses}/{subject}_{ses}_acq-mixed_T1map_raw.nii.gz",
    postprocessed="mri_dataset/derivatives/{subject}/{ses}/{subject}_{ses}_acq-mixed_T1map.nii.gz"
  conda:
    "envs/grip-mixed.yaml"
  params:
    t1_low=100,
    t1_high=10000,
  shell:
    "gmri2fem mri mixed-t1map"
      " --SE {input.se}"
      " --IR {input.ir}"
      " --meta {input.meta}"
      " --output {output.raw}"
      " --postprocessed {output.postprocessed}"
      " --t1_low {params.t1_low}"
      " --t1_high {params.t1_high}"

rule setup_reference_image:
  input:
    "mri_dataset/{subject}/ses-01/anat/{subject}_ses-01_T1w.nii.gz"
  output:
    "mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz"
  shell:
    "cp {input} {output}"


ruleorder: register_pre_contrast_se > register_se_to_reference
rule register_pre_contrast_se:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz",
    moving="mri_dataset/{subject}/ses-01/mixed/{subject}_ses-01_acq-mixed_SE-modulus.nii.gz"
  output:
    "mri_processed_data/{subject}/transforms/{subject}_ses-01_acq-mixed.mat"
  container: "docker://jorgenriseth/greedy:v1.3.0a"
  threads: workflow.cores
  shell:
    "greedy -d 3 -a"
    " -i {input.reference} {input.moving}"
    " -o {output}"
    " -ia-image-centers"
    " -dof 6"
    " -m NMI"
    " -threads {threads}"

rule reslice_pre_contrast_se:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz",
    moving="mri_dataset/{subject}/ses-01/mixed/{subject}_ses-01_acq-mixed_SE-modulus.nii.gz",
    transform="mri_processed_data/{subject}/transforms/{subject}_ses-01_acq-mixed.mat",
  output:
    "mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus_registered.nii.gz",
  container: "docker://jorgenriseth/greedy:v1.3.0a"
  conda:
    "envs/grip-mixed.yaml"
  shell:
    "greedy -d 3"
    " -rf {input.reference}"
    " -ri NN"
    " -rm {input.moving} {output}"
    " -r {input.transform}"
    " -threads {threads}"

rule register_se_to_reference:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus_registered.nii.gz",
    moving="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_SE-modulus.nii.gz"
  output:
    "mri_processed_data/{subject}/transforms/{subject}_{ses}_acq-mixed.mat"
  container: "docker://jorgenriseth/greedy:v1.3.0a"
  threads: workflow.cores
  shell:
    "greedy -d 3 -a"
    " -i {input.reference} {input.moving}"
    " -o {output}"
    " -ia-image-centers"
    " -dof 6"
    " -m NCC 5x5x5"
    " -threads {threads}"

rule reslice_mixed_t1map:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz",
    moving="mri_dataset/derivatives/{subject}/{session}/{subject}_{session}_acq-mixed_T1map_raw.nii.gz",
    transform="mri_processed_data/{subject}/transforms/{subject}_{session}_acq-mixed.mat",
  output:
    "mri_processed_data/{subject}/registered/{subject}_{session}_acq-mixed_T1map_raw_registered.nii.gz"
  container: "docker://jorgenriseth/greedy:v1.3.0a"
  shell:
    "greedy -d 3"
    " -rf {input.reference}"
    " -ri NN"
    " -rm {input.moving} {output}"
    " -r {input.transform}"
    " -threads {threads}"

rule run_segmentation:
  input:
    seg="mri_processed_data/freesurfer/{subject}/mri/aparc+aseg.mgz",
    mask_image="mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus_registered.nii.gz",
  output:
    "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_mask.nii.gz",
    "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_sas-ventricle-cerebellum.nii.gz",
    "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_grouped.nii.gz",
  conda:
    "envs/grip-mixed.yaml"
  shell:
    "python main.py {wildcards.subject}"

rule estimate_csf_concentration:
  input:
    image="mri_processed_data/{subject}/registered/{subject}_{session}_acq-mixed_T1map_raw_registered.nii.gz",
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_T1map_raw_registered.nii.gz",
    mask="mri_processed_data/{subject}/segmentations/{subject}_seg-csf_mask.nii.gz"
  output:
    "mri_processed_data/{subject}/concentrations/{subject}_{session}_acq-mixed_concentration.nii.gz"
  conda:
    "envs/grip-mixed.yaml"
  shell:
    "gmri2fem mri concentration"
    " --input {input.image}"
    " --reference {input.reference}"
    " --output {output}"
    " --mask {input.mask}"
    " --r1 0.004"

rule T1_to_R1:
  input:
    "{anyprefix}T1map{anysuffix}"
  output:
    "{anyprefix}R1map{anysuffix}"
  params:
    t1_low=100,
    t1_high=10000,
  shell:
    "gmri2fem mri t1-to-r1 --input {input} --output {output} --T1_low {params.t1_low} --T1_high {params.t1_high}"
