#!/usr/bin/env bash

base_dir=$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
steinbock="docker run -v ${base_dir}:/data -v /tmp/.X11-unix:/tmp/.X11-unix -v ${HOME}/.Xauthority:/home/steinbock/.Xauthority:ro -e DISPLAY ghcr.io/bodenmillergroup/steinbock:0.8.0"

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
$steinbock segment cellprofiler run

$steinbock measure intensities
$steinbock measure regionprops
$steinbock measure neighbors --type borders --dmax 4
$steinbock measure cellprofiler prepare
$steinbock measure cellprofiler run

$steinbock export ome
$steinbock export csv intensities regionprops -o cells.csv
$steinbock export fcs intensities regionprops -o cells.fcs
$steinbock export anndata --x intensities --obs regionprops
$steinbock export graphs --data intensities --data regionprops
