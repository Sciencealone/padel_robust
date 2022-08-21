#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# v.0.1
# PaDEL-robust: robust version of PaDEL-descriptor Python interface
# Developed in 2022 by Mikhail Markovsky <m.markovsky@gmail.com>


from padel_robust import PadelDescriptor

if __name__ == '__main__':

    # Definition of target substances
    SOURCE_SINGLE = 'CCC'
    SOURCE_ARRAY = [
        'CCC',
        'CCCCC',
        'CCCCCCC'
    ]

    # PaDEL object initiated
    padel = PadelDescriptor(
        convert_3d=True,
        d_2d=True,
        d_3d=True,
        timeout=120,
        threads=1
    )

    # Example of a single calculation
    print('PaDEL-robust test.')
    print('Calculating a single descriptor...')
    descriptors_single = padel.make_descriptor(SOURCE_SINGLE)
    print(f'{len(descriptors_single)} descriptors calculated for a substance.\n')

    # Example of a batch calculation
    print('Starting a batch calculation...')
    descriptors_multiple = padel.make_descriptors_batch(SOURCE_ARRAY)
    print(f'{len(descriptors_multiple[0])} descriptors calculated for each of {len(SOURCE_ARRAY)} substances.\n')
