from symbol import atom_expr
from openmm.openmm import NonbondedForce, CustomBondForce
from openmm.app import ForceField
from openmm.app import forcefield as ff
from openmm.app.forcefield import *

class LJOverrideGenerator(object):
    def __init__(self, forcefield):
        self.ff = forcefield
        self.nbfixTypes = {}
        self.ljTypes = ForceField._AtomTypeParameters(forcefield, 'LJOverrides', 'Atom', ('sigma', 'epsilon'))

    def registerNBFIX(self, parameters):
        types = self.ff._findAtomTypes(parameters, 2)
        if None not in types:
            for type1 in types[0]:
                for type2 in types[1]:
                    epsilon = ff._convertParameterToNumber(parameters['epsilon'])
                    sigma = ff._convertParameterToNumber(parameters['sigma'])
                    self.nbfixTypes[(type1, type2)] = [sigma, epsilon]
                    self.nbfixTypes[(type2, type1)] = [sigma, epsilon]

    def registerLennardJones(self, parameters):
        self.ljTypes.registerAtom(parameters)

    @staticmethod
    def parseElement(element, ff):
        print("IN PARSE")
        existing = [f for f in ff._forces if isinstance(f, LJOverrideGenerator)]
        if len(existing) == 0:
            generator = LJOverrideGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]
        generator.xml = element
        for LJ in element.findall('Atom'):
            generator.registerLennardJones(LJ.attrib)
        for Nbfix in element.findall('NBFixPair'):
            generator.registerNBFIX(Nbfix.attrib)

    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        # First derive the lookup tables.  We need to include entries for every type
        # that a) appears in the system and b) has unique parameters.

        nbfixTypeSet = set().union(*self.nbfixTypes)
        allTypes = set(data.atomType[atom] for atom in data.atoms if data.atomType[atom] in self.ljTypes.paramsForType)
        mergedTypes = []
        mergedTypeParams = []
        paramsToMergedType = {}
        typeToMergedType = {}
        for t in allTypes:
            #if t not in self.ljTypes.paramsForType: continue
            typeParams = self.ljTypes.paramsForType[t]
            params = (typeParams['sigma'], typeParams['epsilon'])
            if t in nbfixTypeSet:
                # NBFIX types cannot be merged.
                typeToMergedType[t] = len(mergedTypes)
                mergedTypes.append(t)
                mergedTypeParams.append(params)
            elif params in paramsToMergedType:
                # We can merge this with another type.
                typeToMergedType[t] = paramsToMergedType[params]
            else:
                # This is a new type.
                typeToMergedType[t] = len(mergedTypes)
                paramsToMergedType[params] = len(mergedTypes)
                mergedTypes.append(t)
                mergedTypeParams.append(params)

        #   default type for those not included
        typeToMergedType['_OTHER_'] = len(mergedTypes)
        mergedTypes.append('_OTHER_')
        mergedTypeParams.append((1.0, 0.0))


        # Now everything is assigned. Create the A- and B-coefficient arrays
        numLjTypes = len(mergedTypes)
        acoef = [0]*(numLjTypes*numLjTypes)
        bcoef = acoef[:]
        for m in range(numLjTypes):
            for n in range(numLjTypes):
                pair = (mergedTypes[m], mergedTypes[n])
                if pair in self.nbfixTypes:
                    sig_old = 0.5*(mergedTypeParams[m][0]+mergedTypeParams[n][0])
                    eps_old = math.sqrt(mergedTypeParams[m][1]*mergedTypeParams[n][1])
                    eps_new = self.nbfixTypes[pair][1]
                    sig_new = self.nbfixTypes[pair][0]
                    A_old, B_old = 4*eps_old*sig_old**12, 4*eps_old*sig_old**6
                    A_new, B_new = 4*eps_new*sig_new**12, 4*eps_new*sig_new**6
                    acoef[m+numLjTypes*n] = A_new - A_old
                    bcoef[m+numLjTypes*n] = B_new - B_old
                    continue
                else:
                    acoef[m+numLjTypes*n] = 0.0
                    bcoef[m+numLjTypes*n] = 0.0

        self.force = mm.CustomNonbondedForce('acoef(type1, type2)/r^12 - bcoef(type1, type2)/r^6;')
        self.force.addTabulatedFunction('acoef', mm.Discrete2DFunction(numLjTypes, numLjTypes, acoef))
        self.force.addTabulatedFunction('bcoef', mm.Discrete2DFunction(numLjTypes, numLjTypes, bcoef))
        self.force.addPerParticleParameter('type')
        self.force.setName('LJOverride')
        if nonbondedMethod in [CutoffPeriodic, Ewald, PME, LJPME]:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        elif nonbondedMethod is NoCutoff:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
        elif nonbondedMethod is CutoffNonPeriodic:
            self.force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffNonPeriodic)
        else:
            raise AssertionError('Unrecognized nonbonded method [%s]' % nonbondedMethod)
        if args['switchDistance'] is not None:
            self.force.setUseSwitchingFunction(True)
            self.force.setSwitchingDistance(args['switchDistance'])
        

        # Add the particles
        # print(acoef)
        # print(bcoef)
        typesToIndicies = dict(zip(mergedTypes, [[] for i in range(numLjTypes)]))
        for atom in data.atoms:
            if data.atomType[atom] in typeToMergedType:
                self.force.addParticle((typeToMergedType[data.atomType[atom]],))
                typesToIndicies[data.atomType[atom]].append(atom.index)
            else:
                self.force.addParticle((typeToMergedType['_OTHER_'],))
            #typesToIndicies[atomType].append(atom.index)

        #   add interaction groups for each NBFix
        addedPairs = []
        for t1, g1 in typesToIndicies.items():
            for t2, g2 in typesToIndicies.items():
                if (t1, t2) in self.nbfixTypes and (t1, t2) not in addedPairs:
                    self.force.addInteractionGroup(g1, g2)
                    addedPairs.append((t1, t2))
                    addedPairs.append((t2, t1))


        self.force.setCutoffDistance(nonbondedCutoff)
        sys.addForce(self.force)

        bondIndices = ff._findBondsForExclusions(data, sys)
        self.force.createExclusionsFromBonds(bondIndices, 3)


        

ff.parsers['LJOverrides'] = LJOverrideGenerator.parseElement