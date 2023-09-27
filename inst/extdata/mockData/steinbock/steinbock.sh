#!/usr/bin/env bash

base_dir=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
steinbock="docker run -v ${base_dir}:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ${HOME}/.Xauthority:/home/steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.16.1"

cd "${base_dir}"

cp -r ../raw .
cp ../panel.csv raw
$steinbock preprocess imc panel
$steinbock preprocess imc images --hpf 50

$steinbock classify ilastik prepare --cropsize 50 --seed 123
cp ../ilastik.ilp pixel_classifier.ilp
rm -r ilastik_crops && mkdir ilastik_crops && cp -r ../analysis/ilastik/*.h5 ilastik_crops
$steinbock classify ilastik fix --no-backup
$steinbock classify ilastik run

$steinbock segment cellprofiler prepare
$steinbock segment cellprofiler run -o masks_ilastik

$steinbock measure intensities --masks masks_ilastik -o intensities_ilastik
$steinbock measure regionprops --masks masks_ilastik -o regionprops_ilastik
$steinbock measure neighbors --type expansion --dmax 4 --masks masks_ilastik -o neighbors_ilastik
$steinbock measure cellprofiler prepare --masks masks_ilastik
$steinbock measure cellprofiler run

# deep learning-based segmentation
$steinbock segment deepcell --minmax -o masks_deepcell

# measurement
$steinbock measure intensities --masks masks_deepcell -o intensities_deepcell
$steinbock measure regionprops --masks masks_deepcell -o regionprops_deepcell
$steinbock measure neighbors --masks masks_deepcell --type expansion --dmax 4 -o neighbors_deepcell

$steinbock export ome
$steinbock export csv intensities_deepcell regionprops_deepcell -o cells.csv
$steinbock export fcs intensities_deepcell regionprops_deepcell -o cells.fcs
$steinbock export anndata --intensities intensities_deepcell --data regionprops_deepcell --neighbors neighbors_deepcell -o cells.h5ad
$steinbock export graphs --data intensities_deepcell --data regionprops_deepcell --neighbors neighbors_deepcell
