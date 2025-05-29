from pathlib import Path
SUBJECTS = sorted(p.name for p in Path("mri_dataset").glob("sub-*"))
rule all:
  input:
    expand(
        "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_grouped.nii.gz",
        subject=SUBJECTS
    )

rule recon_all:
  input:
    t1w="mri_dataset/{subject}/ses-01/sub-01_ses-01_T1w.nii.gz",
    flair="mri_dataset/{subject}/ses-01/sub-01_ses-01_FLAIR.nii.gz",
  output: # protected to avoid accidental reruns.
    protected("mri_processed_data/freesurfer/{{subject}}/mri/aparc+aseg.mgz")
  shell:
    "recon-all -all"
    " -s {wildcards.subject}"
    " -sd mri_processed_data/freesurfer"
    " -i {input.t1w}"
    " -FLAIR {input.flair}"
    " -FLAIRpial"
    " -parallel"

rule estimate_T1maps:
  input:
    se="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_SE-modulus.nii.gz",
    ir="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_IR-corrected-real.nii.gz",
    meta="mri_dataset/{subject}/{ses}/mixed/{subject}_{ses}_acq-mixed_meta.json",
  output:
    "mri_dataset/derivatives/{subject}/{ses}/{subject}_{ses}_acq-mixed_T1map_raw.nii.gz"
  shell:
    "gmri2fem mri t1map-mixed"
      " --SE {input.se}"
      " --IR {input.ir}"
      " --meta {input.meta}"
      " --output {output}"

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
  shell:
    "greedy -d 3"
    " -rf {input.reference}"
    " -ri NN"
    " -rm {input.moving} {output}"
    " -r {input.transform}"

rule register_se_to_reference:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus.nii.gz",
    moving="mri_dataset/{subject}/ses-01/mixed/{subject}_{ses}_acq-mixed_SE-modulus.nii.gz"
  output:
    "mri_processed_data/{subject}/transforms/{subject}_{ses}_acq-mixed.mat"
  shell:
    "greedy -d 3 -a"
    " -i {input.reference} {input.moving}"
    " -o {output}"
    " -ia-image-centers"
    " -dof 6"
    " -m NCC 5x5x5"

rule reslice_mixed_t1map:
  input:
    reference="mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz",
    moving="mri_dataset/derivatives/{subject}/{session}/mixed/{subject}_{session}_acq-mixed_T1map_raw.nii.gz",
    transform="mri_processed_data/{subject}/transforms/{subject}_ses-01_acq-mixed.mat",
  output:
    "mri_processed_data/{subject}/registered/{subject}_{session}_acq-mixed.mat"
  shell:
    "greedy -d 3"
    " -rf {input}"
    " -ri NN"
    " -rm {input.reference} {input.moving}"
    " -r {input.transform}"

rule run_segmentation:
  input:
    seg="mri_processed_data/freesurfer/{subject}/mri/aparc+aseg.mgz",
    mask_image="mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus_registered.nii.gz",
  output:
    "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_sas-ventricle-cerebellum.nii.gz",
    "mri_processed_data/{subject}/segmentations/{subject}_seg-csf_grouped.nii.gz",
  shell:
    "python main.py {wildcards.subject}"

