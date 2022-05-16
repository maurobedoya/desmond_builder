###! /usr/bin/bash
python desmond_builder.py -i config_NVT_ter_langevin_bar_none.dat
python desmond_builder.py -i config_NVT_ter_nose_bar_none.dat
python desmond_builder.py -i config_NPT_ter_langevin_bar_langevin.dat
python desmond_builder.py -i config_NPT_ter_nose_bar_MTK.dat
python desmond_builder.py -i config_NPT_ter_nose_bar_MTK_add_4stages.dat
