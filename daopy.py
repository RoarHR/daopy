import re
import numpy as np
from scipy.spatial import KDTree

N1Names = [f'070{n}' for n in np.arange(184, 194)]
N2Names = [f'080{n}' for n in np.arange(153, 163)]
allNames = N1Names + N2Names

directory = '/home/roar/Documents/Projects/NOT/kic_2569488_again'

doThings = False

def space_str(string, length):
    return ' ' * (length - len(string)) + string

# Find
if doThings:
    with open(f'{directory}/script_find.txt', 'w') as file:
        for name in allNames:
            file.write('at\n')
            file.write(f'{name}.fits\n')
            file.write('fi\n')
            file.write('1 1\n')
            file.write('\n')
            file.write('\n')
            file.write('y\n')

if doThings:
    with open(f'{directory}/script_photometry.txt', 'w') as file:
        for name in allNames:
            file.write('at\n')
            file.write(f'{name}.fits\n')
            file.write('ph\n')
            file.write('\n')
            file.write('\n')
            file.write('\n')
            file.write('\n')

if doThings:
    with open(f'{directory}/script_daomatch.txt', 'w') as file:
        file.write(f'{allNames[0]}.ap\n')
        file.write('\n')
        for name in allNames:
            file.write(f'{name}.ap\n')
        file.write('\n')


def load_ap(file):
    ap = np.loadtxt(file, skiprows=2)
    return np.dstack((ap[0::2, :], ap[1::2, :]))


def reg_coordinates(regfile):
    # Define the regex pattern to find the coordinates
    pattern = re.compile(r'point\((\d+\.?\d*),(\d+\.?\d*)\)')

    x_coords = []
    y_coords = []

    with open(regfile, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                x, y = map(float, match.groups())
                x_coords.append(x)
                y_coords.append(y)

    return np.array(x_coords), np.array(y_coords)

def reg_to_psfstars(regfile, mtrfile):
    # regfile = 'ds9_070184.reg'
    # mtrfile = '070184.mtr'
    mtr = np.loadtxt(f'{directory}/{mtrfile}', skiprows=2)
    mtrTree = KDTree(mtr[:, 1:3])
    x, y = reg_coordinates(f'{directory}/{regfile}')
    _, mtrIndexes = mtrTree.query(np.transpose((x, y)))
    psfStars = mtr[mtrIndexes, 0]
    return psfStars

newPsfStars = reg_to_psfstars('ds9_070184.reg', '070184.mtr')


# psfStars = [1098, 1152, 1123, 780, 1160, 449, 1524, 457]
psfStars = newPsfStars

if doThings:
    for name in allNames:
        mtr = np.loadtxt(f'{directory}/{name}.mtr', skiprows=2)
        ap = load_ap(f'{directory}/{name}.ap')
        mtrRows = np.searchsorted(mtr[:, 0], psfStars)
        apNumbers = mtr[mtrRows, 6].astype(int)
        apRows = np.searchsorted(ap[:, 0, 0], apNumbers)

        with open(f'{directory}/{name}.ap', 'r') as apFile:
            line1 = apFile.readline()
            line2 = apFile.readline()
        line2 = list(line2)
        line2[2] = '3'
        line2 = ''.join(line2)

        with open(f'{directory}/{name}.lst', 'w') as file:
            file.write(line1)
            file.write(line2)
            file.write('\n')
            for row in apRows:
                file.write(space_str(f'{int(ap[row, 0, 0])}', 7))
                file.write(space_str(f'{ap[row, 1, 0]}', 9))
                file.write(space_str(f'{ap[row, 2, 0]}', 9))
                file.write(space_str(f'{ap[row, 3, 0]}', 9))
                file.write(space_str(f'{ap[row, 3, 1]}', 9))
                file.write(space_str(f'0.100\n', 9))

if doThings:
    with open(f'{directory}/script_psf2.txt', 'w') as file:
        for name in allNames:
            file.write('at\n')
            file.write(f'{name}s1.fits\n')
            file.write('ps\n')
            file.write(f'{name}.ap\n')  # File with aperture results
            file.write(f'{name}.lst\n')  # File with PSF stars
            file.write(f'{name}s1.psf\n')  # File for the PSF
            file.write('\n')
            file.write('\n')

if doThings:
    with open(f'{directory}/script_substar2.txt', 'w') as file:
        for name in allNames:
            file.write('at\n')
            file.write(f'{name}.fits\n')
            file.write('su\n')
            file.write(f'{name}s1.psf\n')
            file.write(f'{name}.ap\n')
            file.write('y\n')
            file.write(f'{name}.lst\n')
            file.write(f'{name}s1\n')
            file.write('y\n')

i = 0
if doThings:
    with open(f'{directory}/script_allstar.txt', 'w') as file:
        name = allNames[i]
        file.write(f'\n')
        file.write(f'{name}\n')
        file.write(f'\n')
        file.write(f'\n')
        file.write(f'\n')
        file.write(f'\n')
        i += 1

if doThings:
    with open(f'{directory}/script_new_photometry.txt', 'w') as file:
        for name in allNames:
            file.write('at\n')
            file.write(f'{name}.fits\n')
            file.write('ph\n')
            file.write('\n')
            file.write('\n')
            file.write(f'{name}.coo\n')
            file.write(f'{name}.coo\n')
            file.write(f'\n')
            file.write(f'\n')


