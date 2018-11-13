from __future__ import absolute_import, division, print_function, unicode_literals
import os
import glob
import shutil
import unittest
import numpy as np
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as unit
import mdtraj as md
import AdaptivePELE.adaptiveSampling as adaptiveSampling
from AdaptivePELE.simulation.openmm_simulations import XTCReporter


class TestMD(unittest.TestCase):

    def check_succesful_simulation(self, output, epochs, nTrajs):
        for epoch in range(epochs):
            self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "summary.txt")))
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "trajectory*"))), nTrajs)
            self.assertTrue(len(glob.glob(os.path.join(output, "%d" % epoch, "report*"))), nTrajs)
        self.assertTrue(os.path.exists(os.path.join(output, "%d" % epoch, "clustering", "object.pkl")))

    def testOpenMM3ptb(self):
        output_path = "tests/data/openmm_3ptb"
        controlFile = "tests/data/templetized_controlFile_3ptb_md.conf"

        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testOpenMM1ab1(self):
        output_path = "tests/data/openmm_1ab1"
        controlFile = "tests/data/templetized_controlFile_1ab1_md.conf"
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testRestartAt0(self):
        output_path = "tests/data/openmm_restart_0"
        controlFile = "tests/data/templetized_controlFile_restart_0_md.conf"
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        shutil.copytree("tests/data/restart_0", output_path)
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testRestartAt1(self):
        output_path = "tests/data/openmm_restart_1"
        controlFile = "tests/data/templetized_controlFile_restart_1_md.conf"
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        shutil.copytree("tests/data/restart_1", output_path)
        adaptiveSampling.main(controlFile)
        self.check_succesful_simulation(output_path, 2, 4)
        # cleanup
        shutil.rmtree(output_path)

    def testXTCreporter(self):
        output_PDB = "tests/data/test_xtcreporter.pdb"
        output_XTC = "tests/data/test_xtcreporter.xtc"
        top_PDB = "tests/data/top_xtcreporter.pdb"
        PLATFORM = mm.Platform_getPlatformByName(str('CPU'))
        prmtop = app.AmberPrmtopFile("tests/data/complex.prmtop")
        inpcrd = app.AmberInpcrdFile("tests/data/complex.inpcrd")
        system = prmtop.createSystem(nonbondedMethod=app.PME,
                                     nonbondedCutoff=9*unit.angstroms, constraints=app.HBonds)
        system.addForce(mm.AndersenThermostat(300*unit.kelvin, 1/unit.picosecond))
        integrator = mm.VerletIntegrator(2*unit.femtoseconds)
        force = mm.CustomExternalForce(str("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"))
        force.addGlobalParameter(str("k"), 5.0*unit.kilocalories_per_mole/unit.angstroms**2)
        force.addPerParticleParameter(str("x0"))
        force.addPerParticleParameter(str("y0"))
        force.addPerParticleParameter(str("z0"))
        for j, atom in enumerate(prmtop.topology.atoms()):
            if (atom.name in ('CA', 'C', 'N', 'O') and atom.residue.name != "HOH") or (atom.residue.name == "BEN" and atom.element.symbol != "H"):
                force.addParticle(j, inpcrd.positions[j].value_in_unit(unit.nanometers))
        system.addForce(force)
        simulation = app.Simulation(prmtop.topology, system, integrator, PLATFORM)
        if inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        simulation.context.setPositions(inpcrd.positions)
        simulation.minimizeEnergy(maxIterations=10)
        print("Minimization ended")
        xtcReporter = XTCReporter(output_XTC, 1)
        simulation.reporters.append(app.PDBReporter(output_PDB, 1))
        simulation.reporters.append(xtcReporter)
        simulation.step(10)
        # the XTCReporter does not close the file, so opening the file again
        # without exiting the function causes problems to the mdtraj reader
        xtcReporter.close()
        t_xtc = md.load(str(output_XTC), top=top_PDB)
        t_pdb = md.load(output_PDB)
        self.assertEqual(t_pdb.top, t_xtc.top)
        self.assertEqual(np.sum(np.abs(t_pdb.xyz-t_xtc.xyz) > 1e-3), 0)
        os.remove(output_PDB)
        os.remove(output_XTC)
