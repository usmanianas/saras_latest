#!/usr/bin/python

#############################################################################################################################################
 # Saras
 # 
 # Copyright (C) 2019, Mahendra K. Verma
 #
 # All rights reserved.
 # 
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions are met:
 #     1. Redistributions of source code must retain the above copyright
 #        notice, this list of conditions and the following disclaimer.
 #     2. Redistributions in binary form must reproduce the above copyright
 #        notice, this list of conditions and the following disclaimer in the
 #        documentation and/or other materials provided with the distribution.
 #     3. Neither the name of the copyright holder nor the
 #        names of its contributors may be used to endorse or promote products
 #        derived from this software without specific prior written permission.
 # 
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 # ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 # WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 # DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
 # ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 # (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 # LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 # ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 # SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 #
 ############################################################################################################################################
 ##
 ##! \file checkChannel.py
 #
 #   \brief Python script to compare laminar channel flow results with analytical solution
 #
 #   \author Roshan Samuel
 #   \date Jan 2020
 #   \copyright New BSD License
 #
 ############################################################################################################################################
 ##

ptFile = True

if ptFile:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

from scipy.interpolate import griddata
import numpy as np
import h5py as hp
import yaml as yl

# Pyplot-specific directives
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["mathtext.fontset"] = 'cm'
plt.rcParams["font.weight"] = "medium"

def init():
    global figSize
    global U, V, W
    global Nx, Ny, Nz
    global xLen, yLen, zLen

    Nx = 0
    Ny = 0
    Nz = 0

    xLen = 0.0
    yLen = 0.0
    zLen = 0.0

    U = np.zeros([1, 1, 1])
    V = np.zeros([1, 1, 1])
    W = np.zeros([1, 1, 1])

    figSize = (12, 6)


def parseYAML(paraFile):
    global zBeta
    global Re, dPdX
    global Nx, Ny, Nz
    global xLen, yLen, zLen

    yamlFile = open(paraFile, 'r')
    yamlData = yl.load(yamlFile)

    Nx = 2**yamlData["Mesh"]["X Index"] + 1
    Ny = 2**yamlData["Mesh"]["Y Index"] + 1
    Nz = 2**yamlData["Mesh"]["Z Index"] + 1

    xLen = yamlData["Program"]["X Length"]
    yLen = yamlData["Program"]["Y Length"]
    zLen = yamlData["Program"]["Z Length"]

    Re = yamlData["Program"]["Reynolds Number"]
    dPdX = yamlData["Program"]["Mean Pressure Gradient"]

    zBeta = yamlData["Mesh"]["Z Beta"]


def loadData(timeVal):
    global U, V, W
    global Nx, Ny, Nz

    fileName = "output/Soln_{0:09.4f}.h5".format(float(timeVal))

    try:
        f = hp.File(fileName, 'r')
    except:
        print("Could not open file " + fileName + "\n")
        exit()

    # Initialize and read staggered grid data
    U = np.array(f['Vx'])
    V = np.array(f['Vy'])
    W = np.array(f['Vz'])


def generateZGrid():
    global Nz
    global zLen, zBeta
    global zStag, zColl

    zt = np.linspace(0.0, 1.0, Nz)
    zStag = np.array([zLen*(1.0 - np.tanh(zBeta*(1.0 - 2.0*i))/np.tanh(zBeta))/2.0 for i in zt])
    zColl = np.array([zLen*(1.0 - np.tanh(zBeta*(1.0 - 2.0*i))/np.tanh(zBeta))/2.0 for i in [(zt[j] + zt[j+1])/2 for j in range(len(zt) - 1)]])


def computeAnalytic():
    global zLen
    global zStag
    global Re, dPdX
    global uAnalyticProfile
    global tAnalyticProfile

    delta = zLen/2.0
    uTau = 1.0
    nu = uTau*delta/Re

    rho = 1.0

    # Wall shear stress
    tauw = delta*dPdX

    # Factor for laminar flow profile. This profile is available in many textbooks
    # Please refer Sec 7.1 (Channel Flows) of Turbulence by Pope for more details
    lpFactor = tauw*delta/(2.0*rho*nu)

    # Laminar velocity profile
    normZ = zStag/delta
    uAnalyticProfile = lpFactor*normZ*(2.0 - normZ)

    # Shear stress profile
    tAnalyticProfile = tauw*(1.0 - normZ)


def computeProfiles():
    global Re, zLen
    global U, Nz, zStag
    global uComputedProfile
    global tComputedProfile

    uComputedProfile = np.zeros(Nz)
    for i in range(Nz):
        uComputedProfile[i] = np.mean(U[:, :, i])

    rho = 1.0
    delta = zLen/2.0
    nu = delta/Re
    mu = nu*rho
    tComputedProfile = np.zeros(Nz)
    for i in range(Nz-1):
        tComputedProfile[i] = mu*np.mean((U[:, :, i+1] - U[:, :, i])/(zStag[i+1] - zStag[i]))

    i = Nz-1
    tComputedProfile[i] = mu*np.mean((U[:, :, i] - U[:, :, i-1])/(zStag[i] - zStag[i-1]))


def plotProfiles():
    global U, Nz
    global zStag
    global uAnalyticProfile
    global uComputedProfile
    global tAnalyticProfile
    global tComputedProfile

    # Plot both profiles as subplots in a single figure
    fig, axes = plt.subplots(1, 2, figsize=figSize)

    # Plot mean velocity profile
    axes[0].plot(uAnalyticProfile, zStag, linewidth=2, label='Analytic')
    axes[0].plot(uComputedProfile, zStag, linestyle='--', label='Computed')
    axes[0].set_xlabel(r"$u_x$", fontsize=25)
    axes[0].set_ylabel(r"$z$", fontsize=25)
    axes[0].tick_params(labelsize=20)
    axes[0].legend(fontsize=20)
    axes[0].set_title("Mean Velocity Profile", fontsize=20)

    # Plot shear stress profile
    axes[1].plot(tAnalyticProfile, zStag, linewidth=2, label='Analytic')
    axes[1].plot(tComputedProfile, zStag, linestyle='--', label='Computed')
    axes[1].set_xlabel(r"$\tau$", fontsize=25)
    axes[1].set_ylabel(r"$z$", fontsize=25)
    axes[1].tick_params(labelsize=20)
    axes[1].legend(fontsize=20)
    axes[1].set_title("Shear Stress Profile", fontsize=20)

    plt.gca().set_aspect('auto')
    plt.tight_layout()

    if ptFile:
        plt.savefig("meanProfiles.png")
    else:
        plt.savefig("meanProfiles.png")
        plt.show()


if __name__ == "__main__":
    init()

    parseYAML("input/parameters.yaml")

    loadData(20.0)

    generateZGrid()

    computeAnalytic()
    computeProfiles()

    plotProfiles()

