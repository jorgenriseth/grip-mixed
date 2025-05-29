import os
import click

import simple_mri as sm
import numpy as np
import matplotlib.pyplot as plt

from gmri2fem import segment_tools as segtools
from gmri2fem.csf_segmentation import create_csf_seg
from gmri2fem.visualization import slice_volume
from gmri2fem.segment_tools import find_description_label, find_label_description


@click.command("csf-segmentation")
@click.argument("subject", type=str)
def main(subject):
    aparc_mri = sm.load_mri(
        f"mri_processed_data/freesurfer/{subject}/mri/aparc+aseg.mgz",
        dtype=np.int16,
    )
    mask_base_mri = sm.load_mri(
        f"mri_processed_data/{subject}/registered/{subject}_ses-01_acq-mixed_SE-modulus_registered.nii.gz",
        dtype=np.single,
    )
    csf_seg_mri = create_csf_seg(aparc_mri, mask_base_mri)
    sm.save_mri(
        csf_seg_mri,
        f"mri_processed_data/{subject}/segmentations/{subject}_seg-csf_mask.nii.gz",
        dtype=np.uint8,
    )

    def to_cortex_labels(description_list):
        """Map a cortical region-description to it's corresponding FreeSurfer name
        for both left and right hemisphere."""
        return [f"ctx-lh-{d.lower().replace(' ', '')}" for d in description_list] + [
            f"ctx-rh-{d.lower().replace(' ', '')}" for d in description_list
        ]

    ventricles = [
        "3rd-Ventricle",
        "4th-Ventricle",
        "Left-Lateral-Ventricle",
        "Right-Lateral-Ventricle",
        "Left-Inf-Lat-Vent",
        "Right-Inf-Lat-Vent",
    ]

    infratentorial = [
        "Brain-Stem",
        "Left-Cerebellum-Cortex",
        "Right-Cerebellum-Cortex",
    ]

    basal = [
        *to_cortex_labels(
            [
                "medial orbitofrontal",
                "fusiform",
                "frontal pole",
            ]
        ),
        "Optic-Chiasm",
        "Left-VentralDC",
        "Right-VentralDC",
        "Left-Hippocampus",
        "Right-Hippocampus",
    ]

    medial = to_cortex_labels(
        [
            "Posterior cingulate",
            "Rostral anterior cingulate",
            "Caudal anterior cingulate",
            "Superiorfrontal",
            "Precuneus",
            "Isthmus cingulate",
            "Paracentral",
            "Cuneus",
            "Lingual",
            "Pericalcarine",
            "Entorhinal",
            "Parahippocampal",
        ]
    )

    insula = to_cortex_labels(
        [
            "Insula",
        ]
    )

    groups = {
        "ventricles": ventricles,
        "infratentorial": infratentorial,
        "basal": basal,
        "medial": medial,
        "insula": insula,
    }

    lut_table = segtools.read_lut(
        f"{os.environ['FREESURFER_HOME']}/FreeSurferColorLUT.txt"
    )
    labels = np.unique(csf_seg_mri.data[csf_seg_mri.data > 0])
    grouped_descriptions = sum(groups.values(), start=[])
    grouped_labels = [
        find_description_label(x, lut_table) for x in grouped_descriptions
    ]
    groups["lateral"] = [
        find_label_description(x, lut_table)
        for x in labels[~np.isin(labels, grouped_labels)]
    ]
    relabeling = {
        group: [
            segtools.find_description_label(x, lut_table) for x in group_descriptions
        ]
        for group, group_descriptions in groups.items()
    }

    grouped_seg = segtools.collapse(csf_seg_mri.data, relabeling)
    grouped_seg_mri = sm.SimpleMRI(grouped_seg, csf_seg_mri.affine)
    sm.save_mri(
        grouped_seg_mri,
        f"mri_processed_data/{subject}/segmentations/{subject}_seg-csf_grouped.nii.gz",
        np.int16,
    )

    # Also create a lookup-table for identifying regions
    grouped_lut = segtools.canonical_lut(relabeling, "jet", permute_colors=None)
    grouped_cmap = segtools.listed_colormap(grouped_lut)
    segtools.write_lut(
        f"mri_processed_data/{subject}/segmentations/{subject}_seg-csf_grouped_LUT.txt",
        grouped_lut,
    )

    cerebellum = [
        "Left-Cerebellum-Cortex",
        "Right-Cerebellum-Cortex",
    ]
    labels = np.unique(csf_seg_mri.data[csf_seg_mri.data > 0])
    relabeling = {
        "cerebellum": [find_description_label(x, lut_table) for x in cerebellum],
        "ventricles": [find_description_label(x, lut_table) for x in ventricles],
    }
    relabeling["SAS"] = list(
        labels[~np.isin(labels, relabeling["cerebellum"] + relabeling["ventricles"])]
    )

    groups = {
        desc: [find_label_description(x, lut_table) for x in group_labels]
        for desc, group_labels in relabeling.items()
    }

    newlut = segtools.canonical_lut(relabeling, "viridis", permute_colors=None)
    newseg_cmap = segtools.listed_colormap(newlut)
    newseg = segtools.collapse(csf_seg_mri.data, relabeling)
    newseg_mri = sm.SimpleMRI(newseg, csf_seg_mri.affine)
    sm.save_mri(
        newseg_mri,
        f"mri_processed_data/{subject}/segmentations/{subject}_seg-csf_sas-ventricle-cerebellum.nii.gz",
        np.int16,
    )
    segtools.write_lut(
        f"mri_processed_data/{subject}/segmentations/{subject}_seg-csf_sas-ventricle-cerebellum-LUT.txt",
        newlut,
    )

    lut_table = segtools.read_lut(
        f"{os.environ['FREESURFER_HOME']}/FreeSurferColorLUT.txt"
    )
    fs_cmap = segtools.listed_colormap(lut_table)

    slice_ = ("sagittal", 184)

    t1w_mri = sm.load_mri(
        f"mri_processed_data/{subject}/registered/{subject}_ses-01_T1w_registered.nii.gz",
        dtype=np.single,
    )
    t1w = t1w_mri.data
    im_ref = slice_volume(t1w, *slice_)
    im_segmentations = [
        slice_volume(mri.data, *slice_)
        for mri in [csf_seg_mri, grouped_seg_mri, newseg_mri]
    ]

    fig, axes = plt.subplots(1, 3, figsize=(12, 6))
    for ax in axes:
        ax.imshow(im_ref, cmap="gray", vmin=0, vmax=np.quantile(t1w[t1w > 0], 0.95))
        ax.axis("off")

    for ax, seg, cmap in zip(
        axes, im_segmentations, [fs_cmap, grouped_cmap, newseg_cmap]
    ):
        ax.imshow(np.where(seg > 0, seg, np.nan), **cmap, interpolation="nearest")

    axes[0].set_title("FreeSurfer-based segmentation")
    axes[1].set_title("Grouped segmentation 1")
    axes[2].set_title("Grouped segmentation 2")

    plt.tight_layout()
    plt.savefig("figures/segmentations.png", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
