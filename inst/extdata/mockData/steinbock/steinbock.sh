#!/usr/bin/env bash

base_dir=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
steinbock="docker run -v ${base_dir}:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ${HOME}/.Xauthority:/home/steinbock/.Xauthority:ro -e DISPLAY jwindhager/steinbock:0.3.7"

cd "${base_dir}"

cp -r ../raw .
cp ../panel.csv raw
$steinbock preprocess imc panel
$steinbock preprocess imc images --hpf 50

$steinbock classify ilastik prepare --cropsize 50 --seed 123
cp ../ilastik.ilp pixel_classifier.ilp
rm -r ilastik_crops && cp -r ../analysis/ilastik ilastik_crops
$steinbock classify ilastik fix --no-backup
$steinbock classify ilastik run

$steinbock segment cellprofiler prepare
$steinbock segment cellprofiler run --dest cell_masks

$steinbock measure intensities --mask cell_masks --dest cell_intensities
$steinbock measure regionprops --mask cell_masks --dest cell_regionprops
$steinbock tools data collect cell_intensities --dest cell_intensities.csv
$steinbock tools data collect cell_regionprops --dest cell_regionprops.csv
$steinbock tools data collect cell_intensities cell_regionprops --dest cells.csv

$steinbock measure dists border --mask cell_masks --dest cell_dists
$steinbock measure graphs --dists cell_dists --dmax 4 --dest cell_graphs

$steinbock measure cellprofiler prepare --masks cell_masks
$steinbock measure cellprofiler run

