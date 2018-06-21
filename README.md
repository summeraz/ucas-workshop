# Theoretical and Computational Chemistry: Modeling Methods in Chemical Engineering - Workshop

This repository contains the material required for the workshop portion of the
"Theoretical and Computational Chemistry: Modeling Methods in Chemical Engineering"
course at the University of Chinese Academy of Sciences in June 2018.

There are five workshop sessions designed to onboard students that are new to
molecular simulation to various simulation methods and engines. The MoSDeF toolkit
is utilized as the workhorse for constructing and atom-typing molecular systems.

## Session 1: MoSDeF Tutorial

An overview of the MoSDeF toolkit, focusing primarily on mBuild. This session will
also introduce users to Jupyter notebooks. There are four notebooks associated
with this session:
  - [Introduction to Jupyter notebooks](Session1-MoSDeF_Tutorial/Intro_to_Jupyter.ipynb)
  - [Introduction to MoSDeF](Session1-MoSDeF_Tutorial/Intro_to_MoSDeF.ipynb)
  - [Building simple systems with mBuild](Session1-MoSDeF_Tutorial/Building_an_Alkane.ipynb)
  - [mBuild/Foyer/MD workflow example](Session1-MoSDeF_Tutorial/Workflow_Example.ipynb)

## Session 2: Molecular Dynamics (Lennard-Jones system)

An introduction to molecular dynamics using a simple Lennard-Jones system. This
session will introduce students to the basics of a molecular dynamics simulation
script. The HOOMD engine will be used in this session.

## Session 3: Molecular Dynamics (Atomistic)

An introduction to atomistic molecular dynamics simulation. The GROMACS simulation
engine will be used to simulate systems of bulk liquids.
  - [SPC/E Water](Session3-Atomistic_MD/SPC_E_Water.ipynb)
  - [Bulk alkanes](Session3-Atomistic_MD/Bulk_Alkanes.ipynb)
  - [Alkane/Water mixture](Session3-Atomistic_MD/Alkane_Water_Mixture.ipynb)

## Session 4: Monte Carlo

An introduction to Monte Carlo simulation. Python-wrapped C++ code is used to
perform Monte Carlo simulations of a simple Lennard-Jones system.
  - [Introduction to Monte Carlo](Session4-Monte_Carlo/Monte_Carlo.ipynb)

## Session 5: Non-equilibrium Molecular Dynamics

Advanced molecular dynamics simulations. Specifically, this session will focus on
using the SLLOD equations of motion to perform a non-equilibrium molecular dynamics
simulation of bulk alkanes.
  - [Alkane NEMD](Session5-NEMD/Alkane_Viscosity.ipynb)
