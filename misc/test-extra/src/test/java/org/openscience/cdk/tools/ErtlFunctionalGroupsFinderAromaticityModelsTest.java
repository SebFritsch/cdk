/*
 * To do: Add license header
 */
package org.openscience.cdk.tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.hash.AtomEncoder;
import org.openscience.cdk.hash.BasicAtomEncoder;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.MDLV2000Writer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ErtlFunctionalGroupsFinder.Mode;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;


/**
 * This class can be used to read an SD file containing chemical structures, to extract their functional groups using
 * the ErtlFunctionalGroupsFinder with different settings, and to write the functional groups with their associated 
 * frequency in this SD file to a CSV file.
 * 
 * @author Jonas Schaub
 */
public class ErtlFunctionalGroupsFinderAromaticityModelsTest extends CDKTestCase {
    
    //<editor-fold defaultstate="collapsed" desc="Constants">
    
    //<editor-fold defaultstate="collapsed" desc="Strings">
    /**
     * Called via System.getProperty() to determine a root folder for output directory
     */
    private final String SYSTEM_PROPERTY_FOR_OUTPUT_DIRECTORY = "user.home";
    
    /**
     * Directory for output files
     */
    private final String OUTPUT_DIRECTORY_FROM_SYSTEM_PROPERTY = File.separator + "Documents" + File.separator 
            + "Ertl_fg_finder_test";
    
    /**
     * Resources directory of SD files to be read
     */
    private final String SDF_RESOURCES_PATH = "resources/sdf/";
    
    /**
     * SD file of ChEBI database to be analyzed
     */
    private final String CHEBI_SD_FILE_NAME = "ChEBI_lite.sdf";//ChEBI_lite_3star.sdf";
    
    /**
     * Identifier string for ChEBI test
     */
    private final String CHEBI_TEST_IDENTIFIER = "ChEBI";
    
    /**
     * Identifier string for ChEBI test without counting multiples per molecule
     */
    private final String CHEBI_NO_MULTIPLES_TEST_IDENTIFIER = "ChEBINoMultiples";
    
    /**
     * SD file of ChEMBL database to be analyzed
     */
    private final String CHEMBL_SD_FILE_NAME = "chembl_24.sdf";
    
    /**
     * Identifier string for ChEMBL test
     */
    private final String CHEMBL_TEST_IDENTIFIER = "ChEMBL";
    
    /**
     * Identifier string for test of different CycleFinder settings
     */
    private final String CYCLE_FINDER_TEST_IDENTIFIER = "CycleFinderTest";
    
    /**
     * Name of file for logging occurred exceptions with causing molecules
     */
    private final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log";
    
    /**
     * File type of exceptions log file
     */
    private final String EXCEPTIONS_LOG_FILE_TYPE = ".txt";
    
    /**
     * First line in the exceptions log file
     */
    private final String EXCEPTIONS_LOG_FILE_HEADER = "Following molecules led to the specified exceptions:" 
            + System.getProperty("line.separator")
            + "(Note: If too many exceptions are thrown too fast the JVM stops filling in the complete stack trace. "
            + "You need to be looking at an earlier stack trace to see the details." 
            + System.getProperty("line.separator");
    
    /**
     * Name of file for logging filtered molecules
     */
    private final String FILTERED_MOLECULES_FILE_NAME = "Filtered_Molecules";
    
    /**
     * File type of filtered molecules log file
     */
    private final String FILTERED_MOLECULES_FILE_TYPE = ".txt";
    
    /**
     * First line in the filtered molecules log file
     */
    private final String FILTERED_MOLECULES_FILE_HEADER = "Following molecules were filtered:\n";
    
    /**
     * Format of the time stamp addition to all written output files
     */
    private final String DATE_TIME_FORMAT_PATTERN = "uuuu_MM_dd_HH_mm";
    
    /**
     * Separator for file name segments (test identifier, file name, time stamp)
     */
    private final String FILE_NAME_ADDITION_SEPERATOR = "_";
    
    /**
     * Key for the master HashMap's inner maps under which to store the ChEBI or ChEMBL id of an exemplary molecule 
     * that contains this functional group
     */
    private final String MOLECULE_OF_ORIGIN_KEY = "origin";
    
    //<editor-fold defaultstate="collapsed" desc="Concerning filtering and Preprocessing">
    /**
     * Key for setting an IAtomContainers property that specifies if the molecule consisted of one or more unconnected
     * structures and the biggest of these structures was selected in preprocessing
     */
    private final String BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY = "biggest_fragment_was_selected";
    
    /**
     * Key for setting an IAtomContainers property that specifies if the molecule contained charged atoms and these
     * charges were neutralized in preprocessing
     */
    private final String CHARGES_NEUTRALIZED_PROPERTY_KEY = "charges_were_neutralized";
    
    /**
     * Key for setting an IAtomContainers property that specifies if the molecule must be filtered
     */
    private final String MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY = "molecule_must_be_filtered";
    
    /**
     * Key for setting an IAtomContainers property that specifies why the molecule must be filtered
     */
    private final String CAUSE_FOR_FILTERING_PROPERTY_KEY = "filtering_cause";
    
    /**
     * Message specifying that the atom or bond count of a molecule is zero.
     */
    private final String ATOM_OR_BOND_COUNT_ZERO = "Atom or bond count 0";
    
    /**
     * Message specifying that a molecule contains not allowed atomic numbers.
     */
    private final String FORBIDDEN_ATOMIC_NUMBER = "Contains one or more metal, metalloid or \"R\" atoms";
    //</editor-fold>
    
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder; 
     * String will be split and resulting integers passed to a set
     */
    private final String NON_METALLIC_ATOMIC_NUMBERS = "1,2,6,7,8,9,10,15,16,17,18,34,35,36,53,54,86";
    
    //<editor-fold defaultstate="collapsed" desc="Concerning output file">
    /**
     * Name of file for writing resulting functional groups and frequencies to (output file)
     */
    private final String OUTPUT_FILE_NAME = "Functional_Groups_Frequencies";
    
    /**
     * File type of output file
     */
    private final String OUTPUT_FILE_TYPE = ".csv";
    
    /**
     * Key for the output file's header under which to store the unique SMILES code of a functional group
     */
    private final String SMILES_CODE_KEY = "smiles";
    
    /**
     * Key for the output file's header under which to store the pseudo SMILES code of a functional group
     */
    private final String PSEUDO_SMILES_CODE_KEY = "pseudoSmiles";
    
    /**
     * Key for the output file's header under which to store the hash code code of a functional group
     */
    private final String HASH_CODE_KEY = "hashCode";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdk aromaticity model (and internally for the master HashMap's inner maps)
     */
    private final String CDK_MODEL_SETTINGS_KEY = "cdk";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * daylight aromaticity model (and internally for the master HashMap's inner maps)
     */
    private final String DAYLIGHT_MODEL_SETTINGS_KEY = "daylight";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * piBonds aromaticity model (and internally for the master HashMap's inner maps)
     */
    private final String PIBONDS_MODEL_SETTINGS_KEY = "piBonds";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdkAllowingExocyclic aromaticity model (and internally for the master HashMap's inner maps)
     */
    private final String CDK_EXOCYCLIC_MODEL_SETTINGS_KEY = "cdkExocyclic";
    
    /**
     * Key for the output file's header under which to store the frequency of a functional group when using the 
     * cdk legacy aromaticity model (and internally for the master HashMap's inner maps)
     */
    private final String CDK_LEGACY_MODEL_SETTINGS_KEY = "cdkLegacy";
    
    /**
     * This string will be added to an original settings key when applying the ErtlFunctionalGroupsFinder with activated 
     * generalization
     */
    private final String GENERALIZATION_SETTINGS_KEY_ADDITION = "Generalized";
    
    /**
     * This string will be added to an original settings key (only for exception logging!) when a molecule consists of two or 
     * more unconnected structures and the biggest one is chosen for analysis in the preprocessing
     */
    private final String FRAGMENT_SELECTED_SETTINGS_KEY_ADDITION = "(biggest fragment selected)";
    
    /**
     * This string will be added to an original settings key (only for exception logging!) when a molecule contains charged atoms
     * and these charges are neutralized in the preprocessing
     */
    private final String NEUTRALIZED_SETTINGS_KEY_ADDITION = "(neutralized)";
    
    /**
     * Separator for the output file's values
     */
    private final String OUTPUT_FILE_SEPERATOR = ",";
    
    /**
     * Placeholder String for every functional group's SMILES code whose real SMILES representation could not be 
     * generated
     */
    private final String SMILES_CODE_PLACEHOLDER = "SMILES code could not be created";
    
    /**
     * Placeholder String for 'parent' molecules' ids, assigned in case they could not be read
     */
    private final String MOLECULE_OF_ORIGIN_ID_PLACEHOLDER = "no id for molecule of origin";
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Pseudo-SMILES code">
    /**
     * Pseudo SMILES representation of an aromatic C atom
     */
    private final String PSEUDO_SMILES_AROMATIC_CARBON = "C*";
    
    /**
     * Pseudo SMILES representation of an aromatic N atom
     */
    private final String PSEUDO_SMILES_AROMATIC_NITROGEN = "N*";
    
    /**
     * Pseudo SMILES representation of an aromatic S atom
     */
    private final String PSEUDO_SMILES_AROMATIC_SULPHUR = "S*";
    
    /**
     * Pseudo SMILES representation of an aromatic O atom
     */
    private final String PSEUDO_SMILES_AROMATIC_OXYGEN = "O*";
    
    /**
     * Pseudo SMILES representation of an aromatic Se atom
     */
    private final String PSEUDO_SMILES_AROMATIC_SELENIUM = "Se*";
    
    /**
     * Pseudo SMILES representation of an aromatic P atom
     */
    private final String PSEUDO_SMILES_AROMATIC_PHOSPHOR = "P*";
    
    /**
     * Pseudo SMILES representation of an undefined pseudo atom
     */
    private final String PSEUDO_SMILES_R_ATOM = "R";
    //</editor-fold>
    
    //</editor-fold>
    
    /**
     * Type of the generated SMILES codes
     */
    private final int SMILES_GENERATOR_OUTPUT_MODE = SmiFlavor.Unique;
    
    /**
     * Initial capacity of the master HashMap; May be adjusted when analyzing larger SD files 
     */
    private final int MASTER_HASHMAP_INITIAL_CAPACITY = 10000;
    
    /**
     * Initial capacity of the master HashMap's inner maps that store the frequencies for different settings for a 
     * specific functional group
     */
    private final int INNER_HASHMAPS_INITIAL_CAPACITY = 20;
    
    /**
     * Load factor of the master HashMap
     */
    private final float MASTER_HASHMAP_LOAD_FACTOR = 0.9f;
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Variables">
    /**
     * Directory for all produced files; Will be created by adding SYSTEM_PROPERTY_FOR_OUTPUT_DIRECTORY, 
     * OUTPUT_DIRECTORY_FROM_SYSTEM_PROPERTY and the individual test identifier
     */
    private String outputDirectory;
    
    /**
     * Counts all occurring exceptions in one test
     */
    private int exceptionsCounter;
    
    /**
     * True if the filtered molecules were logged in the filtered molecules log file; This is only necessary in the 
     * first iteration since the applied filters are the same in every iteration (assuming that in a single test
     * only one SD file is analyzed)
     */
    private boolean areFilteredMoleculesLogged;
    
    /**
     * True if all operations in initialize() were successful and the test is able to run
     */
    private boolean isTestAbleToRun;
    
    /**
     * PrintWriter for logging exceptions
     */
    private PrintWriter exceptionsPrintWriter;
    
    /**
     * PrintWriter for logging filtered molecules
     */
    private PrintWriter filteredMoleculesPrintWriter;
    
    /**
     * PrintWriter for writing the resulting functional groups and frequencies
     */
    private PrintWriter dataOutputPrintWriter;
    
    /**
     * SmilesGenerator for generating SMILES and pseudo SMILES representations
     */
    private SmilesGenerator smilesGenerator;
    
    /**
     * MoleculeHashGenerator for easy assessment whether a functional group was already entered into the master HashMap
     */
    private MoleculeHashGenerator molHashGenerator;
    
    /**
     * Instance of the ErtlFunctionalGroupsFinder with generalization turned off
     */
    private ErtlFunctionalGroupsFinder ertlFGFinderGenOff;
    
    /**
     * Instance of the ErtlFunctionalGroupsFinder with generalization turned on
     */
    private ErtlFunctionalGroupsFinder ertlFGFinderGenOn;
    
    /**
     * Master HashMap for storing results; Its keys are the hash codes produced by the MoleculeHashGenerator for the 
     * functional groups and its values are inner HashMaps that hold the IAtomContainer of a functional group and its 
     * frequencies for different settings as String-Object pairs
     */
    private HashMap<Long, HashMap> masterHashMap;
    
    /**
     * A list for storing all used settings keys in a test
     */
    private List<String> settingsKeysList;
    
    /**
     * All allowed atomic numbers to pass to the ErtlFunctionalGroupsFinder as a set of integers (will be parsed from 
     * NON_METALLIC_ATOMIC_NUMBERS)
     */
    private Set<Integer> nonMetallicAtomicNumbersSet;
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * <p>
     * Note: it does not initialize any class variables because that would be unnecessary when it is called by a 
     * test method inherited from CDKTestCase; this is done by initialize()
     */
    public ErtlFunctionalGroupsFinderAromaticityModelsTest() {
        super();
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Tests">
    
    //<editor-fold defaultstate="collapsed" desc="Tests involving databases">
    /**
     * Test for analyzing ChEBI database for all five different aromaticity models supplied by the cdk: daylight, cdk,
     * piBonds, cdkAllowingExocyclic and cdkLegacy
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     *
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occurs
     */
    @Test
    public void analyzeChebi() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CHEBI_TEST_IDENTIFIER, false);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("Loading file named: " + this.CHEBI_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEBI_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 5;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpChebiReader = new IteratingSDFReader(new FileInputStream(tmpChebiSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpChebiReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        //If the 'all' CycleFinder produces an Intractable exception the 'vertexShort' CycleFinder is used
        CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
        
        Aromaticity tmpDaylightModel = new Aromaticity(ElectronDonation.daylight(), tmpCycleFinder);
        Aromaticity tmpCdkModel = new Aromaticity(ElectronDonation.cdk(), tmpCycleFinder);
        Aromaticity tmpPiBondsModel = new Aromaticity(ElectronDonation.piBonds(), tmpCycleFinder);
        Aromaticity tmpCdkAllowingExocyclicModel = new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), tmpCycleFinder);
        Aromaticity tmpCDKLegacyModel = Aromaticity.cdkLegacy();
        
        boolean tmpAreMultiplesCounted = true;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], this.CDK_LEGACY_MODEL_SETTINGS_KEY, tmpCDKLegacyModel, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Test for analyzing ChEBI database for all five different aromaticity models supplied by the cdk: daylight, cdk,
     * piBonds, cdkAllowingExocyclic and cdkLegacy;
     * <p>
     * Difference to analyzeChebi(): If the same functional group occurs multiple times in the same molecule
     * it is counted only once
     *
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occurs
     */
    @Ignore
    @Test
    public void analyzeChebiNoMultiples() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CHEBI_NO_MULTIPLES_TEST_IDENTIFIER, false);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("Loading file named: " + this.CHEBI_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEBI_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 5;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpChebiReader = new IteratingSDFReader(new FileInputStream(tmpChebiSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpChebiReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        //If the 'all' CycleFinder produces an Intractable exception the 'vertexShort' CycleFinder is used
        CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
        
        Aromaticity tmpDaylightModel = new Aromaticity(ElectronDonation.daylight(), tmpCycleFinder);
        Aromaticity tmpCdkModel = new Aromaticity(ElectronDonation.cdk(), tmpCycleFinder);
        Aromaticity tmpPiBondsModel = new Aromaticity(ElectronDonation.piBonds(), tmpCycleFinder);
        Aromaticity tmpCdkAllowingExocyclicModel = new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), tmpCycleFinder);
        Aromaticity tmpCDKLegacyModel = Aromaticity.cdkLegacy();
        
        boolean tmpAreMultiplesCounted = false;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], this.CDK_LEGACY_MODEL_SETTINGS_KEY, tmpCDKLegacyModel, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Test for analyzing ChEMBL database for all five different aromaticity models supplied by the cdk: daylight, cdk,
     * piBonds, cdkAllowingExocyclic and cdkLegacy
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     *
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occurs
     */
    @Ignore
    @Test
    public void analyzeChembl() throws Exception {
        this.initialize(this.CHEMBL_SD_FILE_NAME, this.CHEMBL_TEST_IDENTIFIER, false);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("Loading file named: " + this.CHEMBL_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEMBL_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 5;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpChebiReader = new IteratingSDFReader(new FileInputStream(tmpChebiSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpChebiReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        //If the 'all' CycleFinder produces an Intractable exception the 'vertexShort' CycleFinder is used
        CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
        
        Aromaticity tmpDaylightModel = new Aromaticity(ElectronDonation.daylight(), tmpCycleFinder);
        Aromaticity tmpCdkModel = new Aromaticity(ElectronDonation.cdk(), tmpCycleFinder);
        Aromaticity tmpPiBondsModel = new Aromaticity(ElectronDonation.piBonds(), tmpCycleFinder);
        Aromaticity tmpCdkAllowingExocyclicModel = new Aromaticity(ElectronDonation.cdkAllowingExocyclic(), tmpCycleFinder);
        Aromaticity tmpCDKLegacyModel = Aromaticity.cdkLegacy();
        
        boolean tmpAreMultiplesCounted = true;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], this.CDK_LEGACY_MODEL_SETTINGS_KEY, tmpCDKLegacyModel, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Test for analyzing ChEBI database for six different CycleFinder settings supplied by the cdk: all(), mcb(),
     * relevant(), essential(), tripleShort() and cdkAromaticSet()
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     *
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occurs
     */
    @Ignore
    @Test
    public void analyzeCycleFinderDependency() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CYCLE_FINDER_TEST_IDENTIFIER, false);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("Loading file named: " + this.CHEBI_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEBI_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 6;
        IteratingSDFReader[] tmpReaders = new IteratingSDFReader[tmpRequiredNumberOfReaders];
        try {
            for (int i = 0; i < tmpRequiredNumberOfReaders; i++) {
                IteratingSDFReader tmpChebiReader = new IteratingSDFReader(new FileInputStream(tmpChebiSDFile),
                        DefaultChemObjectBuilder.getInstance(), true);
                tmpReaders[i] = tmpChebiReader;
            }
        } catch (FileNotFoundException aFileNotFoundException) {
            System.out.println("\nSD file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
            return;
        }
        
        Aromaticity tmpDaylightModelAll = new Aromaticity(ElectronDonation.daylight(), Cycles.all());
        Aromaticity tmpDaylightModelMCB = new Aromaticity(ElectronDonation.daylight(), Cycles.mcb());
        Aromaticity tmpDaylightModelRelevant = new Aromaticity(ElectronDonation.daylight(), Cycles.relevant());
        Aromaticity tmpDaylightModelEssential = new Aromaticity(ElectronDonation.daylight(), Cycles.essential());
        Aromaticity tmpDaylightModelTripleShort = new Aromaticity(ElectronDonation.daylight(), Cycles.tripletShort());
        Aromaticity tmpDaylightModelCdkAromaticSet = new Aromaticity(ElectronDonation.daylight(), Cycles.cdkAromaticSet());
        
        boolean tmpAreMultiplesCounted = true;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], "all", tmpDaylightModelAll, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], "mcb", tmpDaylightModelMCB, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], "relevant", tmpDaylightModelRelevant, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], "essential", tmpDaylightModelEssential, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[4], "tripleShort", tmpDaylightModelTripleShort, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[5], "cdkAromaticSet", tmpDaylightModelCdkAromaticSet, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Tests for debugging persisting problems">
    //To do: Find the real reason for clone()'s abnormal behavior!
    /**
     * Test for finding the cause of incorrect valences in the pseudo SMILES code of generalized functional groups.
     *
     * @throws Exception for multiple causes
     */
    /*@Ignore
    @Test
    public void debuggingCloneTest() throws Exception {
        this.initialize("", "DebuggingClone", true);
        String tmpChebi65490Smiles = "Cc1c(O)c2C(=O)C[C@H](Oc2c2C(c3ccccc3)C3=C(Oc12)C(C)(C)C(=O)C(C)(C)C3=O)c1ccccc1";
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        IAtomContainer tmpMolecule = tmpSmilesParser.parseSmiles(tmpChebi65490Smiles);
        MDLV2000Writer tmpMolFileWriter = new MDLV2000Writer(System.out);
        tmpMolecule.removeProperty(CDKConstants.CTAB_SGROUPS);
        tmpMolecule = this.applyFiltersAndPreprocessing(tmpMolecule);
        CycleFinder tmpCycleFinder = Cycles.or(Cycles.all(), Cycles.vertexShort());
        Aromaticity tmpCdkModel = new Aromaticity(ElectronDonation.cdk(), tmpCycleFinder);
        tmpCdkModel.apply(tmpMolecule);
        
        this.ertlFGFinderGenOff.find(tmpMolecule);
        this.ertlFGFinderGenOff.markAtoms(tmpMolecule);
        
        List<IAtomContainer> tmpFunctionalGroups = this.ertlFGFinderGenOff.extractGroups(tmpMolecule);
        for (int i = 0; i < tmpFunctionalGroups.size(); i++) {
            System.out.println();
            System.out.println("----------------------------------------------------------");
            System.out.println();
            System.out.println("Functional group number " + (i+1));
            System.out.println("Pseudo SMILES of functional group: " + this.getPseudoSmilesCode(tmpFunctionalGroups.get(i)));
            System.out.println("Number of implicit Hydrogens: " + AtomContainerManipulator.getImplicitHydrogenCount(tmpFunctionalGroups.get(i)));
            System.out.println("Size of bond array: " + AtomContainerManipulator.getBondArray(tmpFunctionalGroups.get(i)).length);
            System.out.println("Molfile of functional group:");
            tmpMolFileWriter.write(tmpFunctionalGroups.get(i));
            System.out.println();
            List<IAtomContainer> tmpFunctionalGroupsCloned = new LinkedList<>();
            tmpFunctionalGroupsCloned.add(tmpFunctionalGroups.get(i).clone());
            System.out.println("Pseudo SMILES of cloned functional group: " + this.getPseudoSmilesCode(tmpFunctionalGroupsCloned.get(0)));
            System.out.println("Number of implicit Hydrogens: " + AtomContainerManipulator.getImplicitHydrogenCount(tmpFunctionalGroupsCloned.get(0)));
            System.out.println("Size of bond array: " + AtomContainerManipulator.getBondArray(tmpFunctionalGroupsCloned.get(0)).length);
            System.out.println("Molfile of cloned functional group:");
            tmpMolFileWriter.write(tmpFunctionalGroupsCloned.get(0));
            System.out.println();
            List<IAtomContainer> tmpFunctionalGroupsGeneralized = new LinkedList<>();
            tmpFunctionalGroupsGeneralized.add(tmpFunctionalGroups.get(i));
            tmpFunctionalGroupsGeneralized = this.ertlFGFinderGenOff.expandGeneralizedEnvironments(tmpFunctionalGroupsGeneralized);
            System.out.println("Pseudo SMILES of generalized functional group (without prior cloning): " + this.getPseudoSmilesCode(tmpFunctionalGroupsGeneralized.get(0)));
            System.out.println("Number of implicit Hydrogens: " + AtomContainerManipulator.getImplicitHydrogenCount(tmpFunctionalGroupsGeneralized.get(0)));
            System.out.println("Size of bond array: " + AtomContainerManipulator.getBondArray(tmpFunctionalGroupsGeneralized.get(0)).length);
            System.out.println("Molfile of generalized functional group (without prior cloning):");
            tmpMolFileWriter.write(tmpFunctionalGroupsGeneralized.get(0));
            System.out.println();
            List<IAtomContainer> tmpFunctionalGroupsClonedAndGeneralized = this.ertlFGFinderGenOff.expandGeneralizedEnvironments(tmpFunctionalGroupsCloned);
            System.out.println("Pseudo SMILES of generalized functional group with prior cloning: " + this.getPseudoSmilesCode(tmpFunctionalGroupsClonedAndGeneralized.get(0)));
            System.out.println("Number of implicit Hydrogens: " + AtomContainerManipulator.getImplicitHydrogenCount(tmpFunctionalGroupsClonedAndGeneralized.get(0)));
            System.out.println("Size of bond array: " + AtomContainerManipulator.getBondArray(tmpFunctionalGroupsClonedAndGeneralized.get(0)).length);
            System.out.println("Molfile of cloned and then generalized functional group:");
            tmpMolFileWriter.write(tmpFunctionalGroupsClonedAndGeneralized.get(0));
            System.out.println();
            System.out.println();
        }
    }*/
    
    /**
     * Test that shows a resulting functional groups with too many 'R' atoms for input molecules with explicit hydrogens.
     *
     * @throws java.lang.Exception
     */
    @Test public void debuggingTooManyRsTest() throws Exception {
        this.initialize("", "TooManyValences", true);
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        MDLV2000Writer tmpMolFileWriter = new MDLV2000Writer(System.out);
        String tmpAmmoniaSmiles = "N([H])([H])[H]";//CHEBI:16134
        IAtomContainer tmpAmmoniaMol = tmpSmilesParser.parseSmiles(tmpAmmoniaSmiles);
        tmpAmmoniaMol = this.applyFiltersAndPreprocessing(tmpAmmoniaMol);
        Aromaticity.cdkLegacy().apply(tmpAmmoniaMol);
        //System.out.println("Molfile of 'parent' molecule:");
        //tmpMolFileWriter.write(tmpAmmoniaMol);
        List<IAtomContainer> tmpFunctionalGroups = this.ertlFGFinderGenOn.find(tmpAmmoniaMol);
        for (int i = 0; i < tmpFunctionalGroups.size(); i++) {
            System.out.println("----------------------------------------------------------");
            System.out.println();
            System.out.println("Functional group number " + (i+1));
            System.out.println("Pseudo SMILES of functional group: " + this.getPseudoSmilesCode(tmpFunctionalGroups.get(i)));
            System.out.println("Molfile of functional group:");
            tmpMolFileWriter.write(tmpFunctionalGroups.get(i));
            System.out.println("Hash code of functional group: " + this.molHashGenerator.generate(tmpFunctionalGroups.get(i)));
            System.out.println();
        }
        Assert.assertEquals("RN(R)R", this.getPseudoSmilesCode(tmpFunctionalGroups.get(0)));
    }
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Other tests">
    /**
     * Test for correct MoleculeHashGenerator settings/performance on some examples.
     *
     * @throws java.lang.Exception
     */
    @Ignore
    @Test
    public void testMoleculeHashGeneratorSettings() throws Exception {
        this.initialize("", "HashGeneratorTest", true);
        SmilesParser tmpSmilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        
        /*Chebi70986, Chebi16238 and Chebi57692 all contain the same functional group with pseudo SMILES code
        "O=C1N=C(C(=NR)C(=O)N1R)N(R)R", but different hybridizations in the resulting atom containers. But their hash
        codes should be the same under the given settings. This is tested exemplary for many similar cases*/
        String[] tmpSmilesArray = {"OC[C@@H](O)[C@@H](O)[C@@H](O)CN1CC(CO)N=C2C(=O)NC(=O)N=C12",
            "Cc1cc2nc3c(nc(=O)[nH]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C",
            "Cc1cc2nc3c(nc(=O)[n-]c3=O)n(C[C@H](O)[C@H](O)[C@H](O)COP([O-])(=O)OP([O-])(=O)OC[C@H]3O[C@H]([C@H](O)[C@@H]3O)n3cnc4c(N)ncnc34)c2cc1C"};
        List<Long> tmpHashCodesList = new LinkedList<>();
        for (String tmpSmilesCode : tmpSmilesArray) {
            IAtomContainer tmpParsedMolecule = tmpSmilesParser.parseSmiles(tmpSmilesCode);
            tmpParsedMolecule = this.applyFiltersAndPreprocessing(tmpParsedMolecule);
            Aromaticity.cdkLegacy().apply(tmpParsedMolecule);
            List<IAtomContainer> tmpFunctionalGroups = this.ertlFGFinderGenOn.find(tmpParsedMolecule);
            for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroups) {
                if (this.getPseudoSmilesCode(tmpFunctionalGroup).equals("O=C1N=C(C(=NR)C(=O)N1R)N(R)R")) {
                    tmpHashCodesList.add(this.molHashGenerator.generate(tmpFunctionalGroup));
                    /*System.out.println();
                    for (IAtom tmpAtom : tmpFunctionalGroup.atoms()) {
                    System.out.println(tmpAtom.getSymbol() + " : " + tmpAtom.getHybridization());
                    }
                    System.out.println();*/
                }
            }
        }
        for (Long tmpHashCode1 : tmpHashCodesList) {
            for (Long tmpHashCode2 : tmpHashCodesList) {
                Assert.assertEquals(tmpHashCode1.longValue(), tmpHashCode2.longValue());
            }
        }
        
        /*Functional groups like the tertiary amine or the hydroxyl group appear with aromatic and non-aromatic central
        atoms. These two cases should be discrimated by the MoleculeHashGenerator under the given settings*/
        String tmpTertiaryAmineSmiles = "*N(*)*";
        IAtomContainer tmpAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        IAtomContainer tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpTertiaryAmineSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("N"))
                tmpAtom.setIsAromatic(true);
        }
        Assert.assertNotEquals(this.molHashGenerator.generate(tmpAromMol), this.molHashGenerator.generate(tmpNonAromMol));
        String tmpHydroxylGroupSmiles = "[H]O[C]";
        tmpAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        tmpNonAromMol = tmpSmilesParser.parseSmiles(tmpHydroxylGroupSmiles);
        for (IAtom tmpAtom : tmpAromMol.atoms()) {
            if (tmpAtom.getSymbol().equals("C"))
                tmpAtom.setIsAromatic(true);
        }
        Assert.assertNotEquals(this.molHashGenerator.generate(tmpAromMol), this.molHashGenerator.generate(tmpNonAromMol));
        
        /*The following are examples of different (unique!) SMILES codes representing the same functional groups.
        They should be assigned the same hash code*/
        HashMap<String,String> tmpEquivalentSmilesMap = new HashMap<>(20);
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)*", "*N(*)C(=[N](*)*)N(*)*");
        tmpEquivalentSmilesMap.put("*SC1=[N](*)[C]=[C]N1*", "*SC=1N(*)[C]=[C][N]1*");
        tmpEquivalentSmilesMap.put("*[N]1=[C][C]=[C]N1*", "*N1[C]=[C][C]=[N]1*");
        tmpEquivalentSmilesMap.put("*[N](*)=[C]N(*)*", "*N(*)[C]=[N](*)*");
        tmpEquivalentSmilesMap.put("*N(*)[C]=[C][C]=[C][C]=[C][C]=[C][C]=[N](*)*", "*[N](*)=[C][C]=[C][C]=[C][C]=[C][C]=[C]N(*)*");
        tmpEquivalentSmilesMap.put("*[N](*)=C(N(*)*)N(*)P(=O)(O[H])O[H]", "*N(*)C(=[N](*)*)N(*)P(=O)(O[H])O[H]");
        tmpEquivalentSmilesMap.put("[O]I(=O)=O", "O=I(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Br(=O)=O", "O=Br(=O)[O]");
        tmpEquivalentSmilesMap.put("[O]Cl(=O)(=O)=O", "O=Cl(=O)(=O)[O]");
        tmpEquivalentSmilesMap.put("[C]=[C][C]=[C]C#C[C]=[C]C#[C]", "[C]#C[C]=[C]C#C[C]=[C][C]=[C]");
        tmpEquivalentSmilesMap.put("*N1[C]=[C][C]=[N]1*", "*[N]1=[C][C]=[C]N1*");
        for (String tmpKeySmiles : tmpEquivalentSmilesMap.keySet()) {
            IAtomContainer tmpKeyMol = tmpSmilesParser.parseSmiles(tmpKeySmiles);
            IAtomContainer tmpValueMol = tmpSmilesParser.parseSmiles(tmpEquivalentSmilesMap.get(tmpKeySmiles));
            Assert.assertEquals(this.molHashGenerator.generate(tmpKeyMol), this.molHashGenerator.generate(tmpValueMol));
        }
    }
    //</editor-fold>
    
    //</editor-fold>
    
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Initializes all class variables and determines the output directory.
     * 
     * @param aFileName name of the SD file to analyze for a quick pre-check if it is present and the test is therefore 
     * meant to run
     * @param aTestIdentifier a folder with this name will be created in the output directory and it will be added to 
     * the output and log files' names for association of test and files; may be null or empty
     * @param aRunAnyway if true, the test is not aborted even if the file with aFileName can not be found
     * @throws java.lang.Exception if one the FileWriter instances can not be instantiated or more than 
     * Integer.MAX-VALUE tests are to be run this minute (error in the naming of output files) or an unexpected 
     * exception occurs.
     */
    private void initialize(String aFileName, String aTestIdentifier, boolean aRunAnyway) throws Exception {
        System.out.println("\n#########################################################################\n");
        System.out.println("Starting new test, identifier: " + aTestIdentifier);
        System.out.println("\nInitializing class variables...");
        this.isTestAbleToRun = true;
        //First, check if the SD file is present and ignore test if it is not
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile;
        try {
            if (aFileName.equals("") || aFileName == null) {
                aFileName = "[no name was given]";
                throw new NullPointerException();
            }
            tmpSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + aFileName).getFile());
        } catch (NullPointerException aNullPointerException) {
            //getFile() throws a NullPointerException if the required file can not be found
            System.out.println("\n\tUnable to find a file named " + aFileName + " in the resources directory.");
            if (!aRunAnyway) {
                System.out.println("Test is ignored.");
                this.isTestAbleToRun = false;
                Assume.assumeTrue(false);
                return;
            } else {
                System.out.println("\tAs specified, the test will be run anyway.");
            }
        }
        //Determine the output directory
        String tmpOutputRootDirectory = System.getProperty(this.SYSTEM_PROPERTY_FOR_OUTPUT_DIRECTORY);
        if (tmpOutputRootDirectory == null || tmpOutputRootDirectory.isEmpty()) {
            throw new SecurityException("The system property \"" + this.SYSTEM_PROPERTY_FOR_OUTPUT_DIRECTORY + "\" is empty!");
        }
        this.outputDirectory = tmpOutputRootDirectory + this.OUTPUT_DIRECTORY_FROM_SYSTEM_PROPERTY + File.separator + aTestIdentifier;
        File tmpOutputDirectoryFile = new File(this.outputDirectory);
        if (!tmpOutputDirectoryFile.exists()) {
            tmpOutputDirectoryFile.mkdirs();
        }
        System.out.println("\n\tOutput directory: " + this.outputDirectory);
        //Create a time stamp for output and log files
        LocalDateTime tmpDateTime = LocalDateTime.now();
        String tmpDateTimeAddition = tmpDateTime.format(DateTimeFormatter.ofPattern(this.DATE_TIME_FORMAT_PATTERN));
        //Set up exceptions log file
        File tmpExceptionsLogFile = new File(this.outputDirectory + File.separator + aTestIdentifier 
                + this.FILE_NAME_ADDITION_SEPERATOR + this.EXCEPTIONS_LOG_FILE_NAME 
                + this.FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + this.EXCEPTIONS_LOG_FILE_TYPE);
        int tmpFilesInThisMinuteCounter = 1;
        boolean tmpNumberAddedToFileName = false;
        //No pre-existing file is overwritten
        if (tmpExceptionsLogFile.exists()) {
            tmpNumberAddedToFileName = true;
            while (tmpFilesInThisMinuteCounter <= Integer.MAX_VALUE) {
                tmpExceptionsLogFile = new File(this.outputDirectory + File.separator + aTestIdentifier 
                        + this.FILE_NAME_ADDITION_SEPERATOR + this.EXCEPTIONS_LOG_FILE_NAME 
                        + this.FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + "(" + tmpFilesInThisMinuteCounter 
                        + ")" + this.EXCEPTIONS_LOG_FILE_TYPE);
                if (!tmpExceptionsLogFile.exists()) {
                    break;
                }
                if (tmpFilesInThisMinuteCounter == Integer.MAX_VALUE) {
                    throw new Exception("More than [Integer.MAX-VALUE] tests are to be run this minute. "
                            + "This test class is not configurated for that.");
                }
                tmpFilesInThisMinuteCounter++;
            }
        }
        FileWriter tmpExceptionsLogFileWriter = new FileWriter(tmpExceptionsLogFile);
        this.exceptionsPrintWriter = new PrintWriter(tmpExceptionsLogFileWriter);
        this.exceptionsPrintWriter.println(this.EXCEPTIONS_LOG_FILE_HEADER);
        this.exceptionsPrintWriter.flush();
        System.out.println("\tExceptions will be written to: " + tmpExceptionsLogFile.getName());
        //Set up filtered molecules log file
        File tmpFilteredMoleculesFile;
        if (tmpNumberAddedToFileName) {
            tmpFilteredMoleculesFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + this.FILE_NAME_ADDITION_SEPERATOR + this.FILTERED_MOLECULES_FILE_NAME 
                    + this.FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + "(" + tmpFilesInThisMinuteCounter + ")" 
                    + this.FILTERED_MOLECULES_FILE_TYPE);
        } else {
            tmpFilteredMoleculesFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + this.FILE_NAME_ADDITION_SEPERATOR + this.FILTERED_MOLECULES_FILE_NAME 
                    + this.FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + this.FILTERED_MOLECULES_FILE_TYPE);
        }
        FileWriter tmpFilteredMoleculesFileWriter = new FileWriter(tmpFilteredMoleculesFile);
        this.filteredMoleculesPrintWriter = new PrintWriter(tmpFilteredMoleculesFileWriter);
        this.filteredMoleculesPrintWriter.println(this.FILTERED_MOLECULES_FILE_HEADER);
        this.filteredMoleculesPrintWriter.flush();
        System.out.println("\tFiltered molecules will be written to: " + tmpFilteredMoleculesFile.getName());
        //Set up output file
        File tmpOutputFile;
        if (tmpNumberAddedToFileName) {
            tmpOutputFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + this.FILE_NAME_ADDITION_SEPERATOR + this.OUTPUT_FILE_NAME 
                    + FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + "(" + tmpFilesInThisMinuteCounter + ")" 
                    + this.OUTPUT_FILE_TYPE);
        } else {
            tmpOutputFile = new File(this.outputDirectory+ File.separator + aTestIdentifier 
                    + this.FILE_NAME_ADDITION_SEPERATOR+ this.OUTPUT_FILE_NAME 
                    + FILE_NAME_ADDITION_SEPERATOR + tmpDateTimeAddition + this.OUTPUT_FILE_TYPE);
        }
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        this.dataOutputPrintWriter = new PrintWriter(tmpOutputFileWriter);
        System.out.println("\tThe absolute functional groups frequencies will be written to: " + tmpOutputFile.getName());
        
        this.smilesGenerator = new SmilesGenerator(this.SMILES_GENERATOR_OUTPUT_MODE);
        this.molHashGenerator = new HashGeneratorMaker()
                .depth(8)
                .elemental()
                //following line is used instead of .orbital() because the atom hybridizations take more information into 
                //account than the bond order sum but that is not required here
                .encode(BasicAtomEncoder.BOND_ORDER_SUM)
                .encode(CustomAtomEncoder.AROMATICITY) //See enum CustomAtomEncoder below
                .molecular();
        this.ertlFGFinderGenOff = new ErtlFunctionalGroupsFinder(Mode.NO_GENERALIZATION);
        this.ertlFGFinderGenOn = new ErtlFunctionalGroupsFinder(Mode.DEFAULT);
        this.masterHashMap = new HashMap(this.MASTER_HASHMAP_INITIAL_CAPACITY, this.MASTER_HASHMAP_LOAD_FACTOR);
        this.settingsKeysList = new LinkedList<>();
        this.exceptionsCounter = 0;
        this.areFilteredMoleculesLogged = false;
        String[] tmpMetalNumbersStrings = this.NON_METALLIC_ATOMIC_NUMBERS.split(",");
        Integer[] tmpMetalNumbersInt = new Integer[tmpMetalNumbersStrings.length];
        for (int i = 0; i < tmpMetalNumbersStrings.length; i++) {
            tmpMetalNumbersInt[i] = Integer.parseInt(tmpMetalNumbersStrings[i]);
        }
        this.nonMetallicAtomicNumbersSet = new HashSet(Arrays.asList(tmpMetalNumbersInt));
        System.out.println("\n\tDone initializing class variables.\n");
    }
    
    /**
     * Does one iteration of loading molecules from an IteratingSDFReader, applying the given aromaticity model, 
     * extracting their functional groups (with and without generalization) and adding the results to the master HashMap.
     * Exceptions caused by read-in molecules are caught and logged.
     * 
     * @param aReader to load the molecules to be screened
     * @param aSettingsKey resulting functional groups will be added to the inner maps of the master HasMap under this key
     * @param anAromaticity to apply to the molecules
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in the same molecule will 
     * only be counted once
     */
    private void calculateAbsoluteFGFrequencies(
            IteratingSDFReader aReader, 
            String aSettingsKey, 
            Aromaticity anAromaticity, 
            boolean anAreMultiplesCounted) {
        
        System.out.println("\nAnalyzing database using specified settings: " + aSettingsKey);
        //<editor-fold defaultstate="collapsed" desc="Counter definitions">
        int tmpMoleculesCounter = 0; //total number of analyzed molecules
        int tmpChargeCounter = 0; //number of molecules that contain one or more charged atoms
        int tmpUnconnectedCounter = 0; //number of structures with unconnected substructures
        int tmpFilteredMoleculesCounter = 0; //filtered molecules, see below
        int tmpMetallicCounter = 0; //number of molecules that were filtered because they contain one or more metal, metalloid or "R" atoms
        int tmpNoAtomOrBondCounter = 0; //number of molecules that were filtered because they do not contain an atom or a bond
        int tmpSkippedMoleculesCounter = 0; //molecules that caused exceptions and are therefore skipped
        int tmpNoFunctionalGroupsCounter = 0; //number of molecules that contain no functional group
        //</editor-fold>
        while (aReader.hasNext()) {
            tmpMoleculesCounter++;
            List<IAtomContainer> tmpFunctionalGroups;
            List<IAtomContainer> tmpFunctionalGroupsGeneralized;
            IAtomContainer tmpMolecule = (IAtomContainer) aReader.next();
            if (tmpMolecule == null) {
                break;
            }
            String tmpSettingsKeyForLogging = aSettingsKey;
            IAtomContainer tmpOriginalMolecule = null;
            try {
                /*Note: A molecule can be supposed to be filtered for more than one of the named reasons but only the first 
                  tested reason will be named as cause*/
                String tmpCauseForFiltering = "";
                //Remove s-groups ... they cause trouble in atomContainer.clone()(->SgroupManilupator.copy)
                tmpMolecule.removeProperty(CDKConstants.CTAB_SGROUPS);
                //Produce a clone of the original molecule for logging
                tmpOriginalMolecule = tmpMolecule.clone();
                //Apply filters and preprocessing
                tmpMolecule = this.applyFiltersAndPreprocessing(tmpMolecule);
                //Filter molecule if necessary
                if (tmpMolecule.getProperty(this.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY)) {
                    tmpCauseForFiltering = tmpMolecule.getProperty(this.CAUSE_FOR_FILTERING_PROPERTY_KEY);
                    if (tmpCauseForFiltering.equals(this.ATOM_OR_BOND_COUNT_ZERO)) {
                        tmpNoAtomOrBondCounter++;
                    } else if (tmpCauseForFiltering.equals(this.FORBIDDEN_ATOMIC_NUMBER)) {
                        tmpMetallicCounter++;
                    }
                    tmpFilteredMoleculesCounter++;
                    if (!this.areFilteredMoleculesLogged) {
                        this.logFilteredMolecule(tmpMolecule, tmpFilteredMoleculesCounter, tmpCauseForFiltering);
                    }
                    continue;
                }
                //Add according additions to the settings key for logging if the molecule was preprocessed
                if (tmpMolecule.getProperty(this.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY)) {
                    tmpUnconnectedCounter++;
                    tmpSettingsKeyForLogging += this.FRAGMENT_SELECTED_SETTINGS_KEY_ADDITION;
                }
                if (tmpMolecule.getProperty(this.CHARGES_NEUTRALIZED_PROPERTY_KEY)) {
                    tmpChargeCounter++;
                    tmpSettingsKeyForLogging += this.NEUTRALIZED_SETTINGS_KEY_ADDITION;
                }
                
                anAromaticity.apply(tmpMolecule);
                //Now the analysis of functional groups
                //throws UnsupportedOperationException!
                //tmpFunctionalGroups = this.ertlFGFinderGenOff.find(tmpMolecule);
                tmpFunctionalGroups = new LinkedList<>(); //Workaround!
                //Do the extraction again with activated generalization
                tmpSettingsKeyForLogging += this.GENERALIZATION_SETTINGS_KEY_ADDITION;
                //clone the original molecule before that or rather before first marking/extraction?
                tmpFunctionalGroupsGeneralized = this.ertlFGFinderGenOn.find(tmpMolecule);
            } catch (Exception anException) {
                tmpSkippedMoleculesCounter++;
                if (tmpOriginalMolecule != null) {
                    this.logException(anException, tmpSettingsKeyForLogging, tmpOriginalMolecule);
                } else {
                    this.logException(anException, tmpSettingsKeyForLogging, tmpMolecule);
                }
                continue;
            }
            //if (!(tmpFunctionalGroups == null || tmpFunctionalGroups.isEmpty())) {
            if (!(tmpFunctionalGroupsGeneralized == null || tmpFunctionalGroupsGeneralized.isEmpty())) {
                //If a molecule does not have FGs without generalization it does not have FGs with generalization either
                this.enterFunctionalGroupsIntoMasterMap(tmpFunctionalGroups, 
                        aSettingsKey, 
                        anAreMultiplesCounted, 
                        tmpOriginalMolecule);
                this.enterFunctionalGroupsIntoMasterMap(tmpFunctionalGroupsGeneralized, 
                        aSettingsKey + this.GENERALIZATION_SETTINGS_KEY_ADDITION,
                        anAreMultiplesCounted, 
                        tmpOriginalMolecule);
            } else {
                tmpNoFunctionalGroupsCounter++;
            }
        }
        try {
            aReader.close();
        } catch (IOException anIOException) { }
        //Since the filters remain the same in every iteration filtered molecules must be logged only once
        //(assuming that only one SD file is analyzed in a test)
        if (!this.areFilteredMoleculesLogged) {
            this.areFilteredMoleculesLogged = true;
        }
        System.out.println("Analysis done");
        System.out.println(tmpMoleculesCounter + " molcules were screened.");
        System.out.println(tmpChargeCounter + " molecules contained one or more charged atom.");
        System.out.println(tmpUnconnectedCounter + " molecules had two or more unconnected structures.");
        System.out.println(tmpFilteredMoleculesCounter + " molecules were filtered.");
        System.out.println("\t" + tmpNoAtomOrBondCounter + " molecules were filtered because their atom or bond count is 0.");
        System.out.println("\t" + tmpMetallicCounter + " molecules were filtered because they contained one or more "
                + "metal, metalloid or \"R\" atom.");
        System.out.println(tmpSkippedMoleculesCounter + " molecules were skipped due to exceptions.");
        System.out.println(tmpNoFunctionalGroupsCounter + " molecules contained no functional groups.");
        System.out.println(this.masterHashMap.size() + " different functional groups were detected so far.");
        
    }
    
    /**
     * Combines all filtering and preprocessing steps. Molecules will be filtered if they contain a not allowed atomic 
     * number (metal, metalloid or 'R' atoms) or if their atom or bond count is zero. If the molecule should be filtered
     * the IAtomContainer property MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY will be set to true and the 
     * CAUSE_FOR_FILTERING_PROPERTY_KEY will give the cause for the filtering.
     * <p>
     * The preprocessing consists of neutralizing any charges in the molecule and selecting the biggest fragment for
     * further processing if the molecule consists of one or more unconnected structures. If any of these cases apply
     * the IAtomContainer properties CHARGES_NEUTRALIZED_PROPERTY_KEY and BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY will be
     * set accordingly.
     * 
     * @param aMolecule the molecule to be processed
     * @return the processed molecule
     * @throws CDKException if AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms() or 
     * applyNeutralizeChargesPreprocessing() throws a CDKException
     */
    private IAtomContainer applyFiltersAndPreprocessing(IAtomContainer aMolecule) throws CDKException {
        aMolecule.setProperty(this.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, false);
        aMolecule.setProperty(this.CAUSE_FOR_FILTERING_PROPERTY_KEY, "");
        aMolecule.setProperty(this.CHARGES_NEUTRALIZED_PROPERTY_KEY, false);
        aMolecule.setProperty(this.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY, false);
        if (this.isAtomOrBondCountZero(aMolecule)) {
            aMolecule.setProperty(this.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, true);
            aMolecule.setProperty(this.CAUSE_FOR_FILTERING_PROPERTY_KEY, this.ATOM_OR_BOND_COUNT_ZERO);
            return aMolecule;
        }
        aMolecule.removeProperty(CDKConstants.CTAB_SGROUPS); //In case it was not done earlier
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(aMolecule);
        //##############################################################################################################
        //Workaround for solving problems with explicit hydrogens in the input molecules
        //AtomContainerManipulator.suppressHydrogens(aMolecule);
        //CDKHydrogenAdder.getInstance(aMolecule.getBuilder()).addImplicitHydrogens(aMolecule);
        //##############################################################################################################
        //Preprocessing: From structures containing two or more unconnected structures (e.g. ions) 
        //choose the largest structure for analysis
        //This step must be done prior to the next step!
        if (this.isStructureUnconnected(aMolecule)) {
            aMolecule = this.applySelectBiggestFragmentPreprocessing(aMolecule);
            aMolecule.setProperty(this.BIGGEST_FRAGMENT_SELECTED_PROPERTY_KEY, true);
        }
        //Filter molecules containing metals, metalloids or "R" atoms
        if (this.areMetallicOrMetalloidAtomsInMolecule(aMolecule)) {
            aMolecule.setProperty(this.MOLECULE_MUST_BE_FILTERED_PROPERTY_KEY, true);
            aMolecule.setProperty(this.CAUSE_FOR_FILTERING_PROPERTY_KEY, this.FORBIDDEN_ATOMIC_NUMBER);
            return aMolecule;
        }
        //Preprocessing: Neutralize charges
        if (this.isMoleculeCharged(aMolecule)) {
            aMolecule = this.applyNeutralizeChargesPreprocessing(aMolecule);
            aMolecule.setProperty(this.CHARGES_NEUTRALIZED_PROPERTY_KEY, true);
        }
        return aMolecule;
    }
    
    /**
     * Returns true, if the atom or bond count of the molecule is zero. This is a cause for filtering the molecule.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the atom or bond count of the molecule is zero
     */
    private boolean isAtomOrBondCountZero(IAtomContainer aMolecule) {
        return aMolecule.getAtomCount() == 0 || aMolecule.getBondCount() == 0;
    }
    
    /**
     * Returns true, if a not allowed atomic number is detected in the molecule. This is a cause for filtering the molecule.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule contains a not allowed element
     */
    private boolean areMetallicOrMetalloidAtomsInMolecule(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (!this.nonMetallicAtomicNumbersSet.contains(tmpAtom.getAtomicNumber())) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Returns true, if the molecule consists of two or more unconnected structures.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule consists of two or more unconnected structures
     */
    private boolean isStructureUnconnected(IAtomContainer aMolecule) {
        return (!ConnectivityChecker.isConnected(aMolecule));
    }
    
    /**
     * Returns the biggest unconnected fragment/structure of the given molecule. To pre-check if the molecule has 
     * two or more unconnected structures use isStructureConnected(). All set properties of aMolecule will be copied to 
     * the returned IAtomContainer.
     * 
     * @param aMolecule the molecule whose biggest fragment should be found
     * @return the biggest unconnected fragment/structure of the given molecule
     */
    private IAtomContainer applySelectBiggestFragmentPreprocessing(IAtomContainer aMolecule) {
        IAtomContainerSet tmpFragmentsSet = ConnectivityChecker.partitionIntoMolecules(aMolecule);
        IAtomContainer tmpBiggestFragment = null;
        for (IAtomContainer tmpFragment : tmpFragmentsSet.atomContainers()) {
            if (tmpBiggestFragment == null || tmpBiggestFragment.getAtomCount() < tmpFragment.getAtomCount()) {
                tmpBiggestFragment = tmpFragment;
            }
        }
        tmpBiggestFragment.setProperties(aMolecule.getProperties());
        return tmpBiggestFragment;
    }
    
    /**
     * Returns true, if the molecule contains charged atoms.
     * 
     * @param aMolecule the molecule to be tested
     * @return true, if the molecule contains charged atoms
     */
    private boolean isMoleculeCharged(IAtomContainer aMolecule) {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getFormalCharge() != 0) {
                return true;
            }
        }
        return false;
    }
    
    /**
     * Neutralizes all non-zero charges in the given molecule. To pre-check if the molecule has charged atoms use 
     * isMoleculeCharged().
     * 
     * @param aMolecule the molecule to be neutralized
     * @return the same IAtomContainer instance as aMolecule but with neutralized charges
     * @throws CDKException if CDKAtomTypeMatcher.findMatchingAtomType() or CDKHydrogenAdder.addImplicitHydrogens 
     * throws a CDKException
     */
    private IAtomContainer applyNeutralizeChargesPreprocessing(IAtomContainer aMolecule) throws CDKException {
        for (IAtom tmpAtom : aMolecule.atoms()) {
            if (tmpAtom.getFormalCharge() != 0) {
                tmpAtom.setFormalCharge(0);
                CDKHydrogenAdder tmpHAdder = CDKHydrogenAdder.getInstance(aMolecule.getBuilder());
                CDKAtomTypeMatcher tmpMatcher = CDKAtomTypeMatcher.getInstance(aMolecule.getBuilder());
                IAtomType tmpMatchedType = tmpMatcher.findMatchingAtomType(aMolecule, tmpAtom);
                if (tmpMatchedType != null) {
                    AtomTypeManipulator.configure(tmpAtom, tmpMatchedType);
                }
                tmpHAdder.addImplicitHydrogens(aMolecule, tmpAtom);
            }
        }
        return aMolecule;
    }
    
    /**
     * Inserts a list of IAtomContainers (the functional groups of one molecule) into the master HashMap. If the 
     * functional group is already inserted (with this settings key) its frequency for the given settings key is raised 
     * by one or else a new inner HashMap will be created for it.
     * 
     * @param aFunctionalGroupsList the functional groups of one molecule to be inserted
     * @param aSettingsKey will be the key of the inner HashMap inside the master HashMap
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in aFunctionalGroupsList will 
     * only be entered once into the master Hashmap
     * @param anFGContainingMolecule the molecule from which the functional groups originated; will be added to the 
     * master Hashmap
     */
    private void enterFunctionalGroupsIntoMasterMap(
            List<IAtomContainer> aFunctionalGroupsList, 
            String aSettingsKey, 
            boolean anAreMultiplesCounted,
            IAtomContainer anFGContainingMolecule) {
        
        if (!this.settingsKeysList.contains(aSettingsKey)) {
            this.settingsKeysList.add(aSettingsKey);
        }
        List<Long> tmpAlreadyEnteredFGsForThisMol = new LinkedList<>();
        for (IAtomContainer tmpFunctionalGroup : aFunctionalGroupsList) {
            long tmpHashCode = this.molHashGenerator.generate(tmpFunctionalGroup);
            //Case: Multiples are counted only once and the functional group was already entered for this molecule
            if (!anAreMultiplesCounted && tmpAlreadyEnteredFGsForThisMol.contains(tmpHashCode)) {
                continue;
            }
            //Case: functional group is already in the master HashMap
            if (this.masterHashMap.containsKey(tmpHashCode)) {
                HashMap<String, Object> tmpInnerMap = (HashMap)this.masterHashMap.get(tmpHashCode);
                //And a key-value pair for this settings key is already present too -> raise frequency by one
                if (tmpInnerMap.containsKey(aSettingsKey)) {
                    int tmpFrequency = (int)tmpInnerMap.get(aSettingsKey);
                    tmpFrequency++;
                    tmpInnerMap.put(aSettingsKey, tmpFrequency);
                //there is no key-value pair for this settings key in the inner HashMap -> create one
                } else {
                    tmpInnerMap.put(aSettingsKey, 1);
                }
            //The functional group did not occur before -> create a new inner HashMap for this molecule
            } else {
                HashMap<String,Object> tmpNewInnerMap = new HashMap(this.INNER_HASHMAPS_INITIAL_CAPACITY);
                tmpNewInnerMap.put(this.MOLECULE_OF_ORIGIN_KEY, anFGContainingMolecule);
                tmpNewInnerMap.put(aSettingsKey, 1);
                String tmpSmilesCode;
                String tmpPseudoSmilesCode;
                try {
                    //Creation of unique SMILES code
                    tmpSmilesCode = this.smilesGenerator.create(tmpFunctionalGroup);
                    //Creation of pseudo SMILES code
                    tmpPseudoSmilesCode = this.getPseudoSmilesCode(tmpFunctionalGroup);
                } catch (CDKException | NullPointerException | CloneNotSupportedException anException) {
                    this.logException(anException, aSettingsKey + "Creating SMILES code", tmpFunctionalGroup);
                    tmpSmilesCode = this.SMILES_CODE_PLACEHOLDER;
                    tmpPseudoSmilesCode = this.SMILES_CODE_PLACEHOLDER;
                }
                tmpNewInnerMap.put(this.SMILES_CODE_KEY, tmpSmilesCode);
                tmpNewInnerMap.put(this.PSEUDO_SMILES_CODE_KEY, tmpPseudoSmilesCode);
                this.masterHashMap.put(tmpHashCode, tmpNewInnerMap);
            }
            tmpAlreadyEnteredFGsForThisMol.add(tmpHashCode);
        }   
    }
    
    /**
     * Writes all frequency data with the respective hash code, SMILES code, pseudo SMILES code and the ChEBI or ChEMBL 
     * id of an exemplary molecule that contains this functional group for all functional groups in the master HashMap 
     * to the output file.
     * <p>
     * Note: The IAtomContainer object stored in the master HashMap's inner maps are cloned in this method for pseudo 
     * SMILES creation. And all PrintWriter instances will be closed.
     */
    private void postProcessAndSaveData() {
        System.out.println("\nWriting to file...");
        //Writing the output file's header
        String tmpFileHeader = this.HASH_CODE_KEY + this.OUTPUT_FILE_SEPERATOR + this.PSEUDO_SMILES_CODE_KEY 
                + this.OUTPUT_FILE_SEPERATOR + this.SMILES_CODE_KEY;
        for (String tmpSettingsKey : this.settingsKeysList) {
            tmpFileHeader += this.OUTPUT_FILE_SEPERATOR + tmpSettingsKey;
        }
        tmpFileHeader += this.OUTPUT_FILE_SEPERATOR + this.MOLECULE_OF_ORIGIN_KEY;
        this.dataOutputPrintWriter.println(tmpFileHeader);
        this.dataOutputPrintWriter.flush();
        Iterator tmpFunctionalGroupsIterator = this.masterHashMap.keySet().iterator();
        //Iteration for all molecules in the master HashMap
        while (tmpFunctionalGroupsIterator.hasNext()) {
            long tmpHashCode = (long)tmpFunctionalGroupsIterator.next();
            HashMap tmpInnerMap = (HashMap)this.masterHashMap.get(tmpHashCode);
            String tmpSmilesCode = (String) tmpInnerMap.get(this.SMILES_CODE_KEY);
            String tmpPseudoSmilesCode = (String) tmpInnerMap.get(this.PSEUDO_SMILES_CODE_KEY);
            //Writing the record for this functional group
            String tmpRecord = tmpHashCode + this.OUTPUT_FILE_SEPERATOR + tmpPseudoSmilesCode 
                    + this.OUTPUT_FILE_SEPERATOR + tmpSmilesCode;
            for (String tmpSettingsKey : this.settingsKeysList) {
                if (tmpInnerMap.get(tmpSettingsKey) == null) {
                    tmpInnerMap.put(tmpSettingsKey, 0);
                }
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + tmpInnerMap.get(tmpSettingsKey);
            }
            IAtomContainer tmpMoleculeOfOrigin = (IAtomContainer)tmpInnerMap.get(this.MOLECULE_OF_ORIGIN_KEY);
            String tmpChebiId = tmpMoleculeOfOrigin.getProperty("ChEBI ID");
            String tmpChemblId = tmpMoleculeOfOrigin.getProperty("chembl_id");
            String tmpCdkTitle = tmpMoleculeOfOrigin.getProperty(CDKConstants.TITLE);
            if (tmpChebiId != null) {
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + tmpChebiId;
            } else if (tmpChemblId != null) {
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + tmpChemblId;
            } else if (tmpCdkTitle != null) {
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + tmpCdkTitle;
            } else {
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + this.MOLECULE_OF_ORIGIN_ID_PLACEHOLDER;
            }
            this.dataOutputPrintWriter.println(tmpRecord);
            this.dataOutputPrintWriter.flush();
            tmpFunctionalGroupsIterator.remove();
        }
        this.dataOutputPrintWriter.close();
        this.exceptionsPrintWriter.close();
        this.filteredMoleculesPrintWriter.close();
    }
    
    /**
     * Returns the Pseudo SMILES code of the given molecule. Pseudo atoms are represented by 'R' atoms and aromatic 
     * (C, S, O, P, Se, N) atoms will be marked by *.
     * 
     * @param aMolecule the molecule to be represented by the Pseudo SMILES string
     * @return the Pseudo SMILES representation of the given molecule
     * @throws CDKException if SmilesGenerator.create() throws a CDKException
     * @throws CloneNotSupportedException if IAtomContainer.clone() throws a CloneNotSupportedException when 
     * invoked on aMolecule
     */
    private String getPseudoSmilesCode(IAtomContainer aMolecule) throws CDKException, CloneNotSupportedException {
        IAtomContainer tmpMolecule = aMolecule.clone();
        for (IAtom tmpAtom : tmpMolecule.atoms()) {
            if (tmpAtom.isAromatic()) {
                IAtom tmpReplacementAtom = null;
                if (null != tmpAtom.getSymbol()) switch (tmpAtom.getSymbol()) {
                    case "C":
                        tmpReplacementAtom = new Atom("Ce");
                        break;
                    case "N":
                        tmpReplacementAtom = new Atom("Nd");
                        break;
                    case "S":
                        tmpReplacementAtom = new Atom("Sm");
                        break;
                    case "O":
                        tmpReplacementAtom = new Atom("Os");
                        break;
                    case "Se":
                        tmpReplacementAtom = new Atom("Sc");
                        break;
                    case "P":
                        tmpReplacementAtom = new Atom("Pm");
                        break;
                    default:
                        break;
                }
                if (tmpReplacementAtom != null) {
                    Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                    AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                    tmpReplacementAtom.setImplicitHydrogenCount(
                            tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                }
            }
            tmpAtom.setIsAromatic(false);
            if (tmpAtom instanceof IPseudoAtom && "R".equals(((IPseudoAtom)tmpAtom).getLabel())) {  
                //second condition: see creation of R atoms in ErtlFunctionalGroupsFinder
                IAtom tmpReplacementAtom = new Atom("Es");
                Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                AtomContainerManipulator.replaceAtomByAtom(tmpMolecule, tmpAtom, tmpReplacementAtom);
                tmpReplacementAtom.setImplicitHydrogenCount(
                        tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
            }
        }
        for (IBond tmpBond : tmpMolecule.bonds()) {
            tmpBond.setIsAromatic(false);
        }
        String tmpPseudoSmilesCode = this.smilesGenerator.create(tmpMolecule);
        tmpPseudoSmilesCode = tmpPseudoSmilesCode
                .replaceAll("(\\[Es\\])", this.PSEUDO_SMILES_R_ATOM)
                .replaceAll("(Es)", this.PSEUDO_SMILES_R_ATOM)
                .replaceAll("(\\[Ce\\])", this.PSEUDO_SMILES_AROMATIC_CARBON)
                .replaceAll("(Ce)", this.PSEUDO_SMILES_AROMATIC_CARBON)
                .replaceAll("(\\[Nd\\])", this.PSEUDO_SMILES_AROMATIC_NITROGEN)
                .replaceAll("(Nd)", this.PSEUDO_SMILES_AROMATIC_NITROGEN)
                .replaceAll("(\\[Sm\\])", this.PSEUDO_SMILES_AROMATIC_SULPHUR)
                .replaceAll("(Sm)", this.PSEUDO_SMILES_AROMATIC_SULPHUR)
                .replaceAll("(\\[Os\\])", this.PSEUDO_SMILES_AROMATIC_OXYGEN)
                .replaceAll("(Os)", this.PSEUDO_SMILES_AROMATIC_OXYGEN)
                .replaceAll("(\\[Sc\\])", this.PSEUDO_SMILES_AROMATIC_SELENIUM)
                .replaceAll("(Sc)", this.PSEUDO_SMILES_AROMATIC_SELENIUM)
                .replaceAll("(\\[Pm\\])", this.PSEUDO_SMILES_AROMATIC_PHOSPHOR)
                .replaceAll("(Pm)", this.PSEUDO_SMILES_AROMATIC_PHOSPHOR);
        return tmpPseudoSmilesCode;
    }
    
    /**
     * Logs molecules that are filtered from the SD file to the filtered molecules file with SMILES code, ChEBI name 
     * and id or ChEMBL id and why they were filtered.
     * 
     * @param aMolecule the filtered molecule to be logged
     * @param aCounter the number of filtered molecules so far (will be written to file)
     * @param aCause why this molecule was filtered
     */
    private void logFilteredMolecule(IAtomContainer aMolecule, int aCounter, String aCause) {
        this.filteredMoleculesPrintWriter.println();
        this.filteredMoleculesPrintWriter.println(aCounter + ". Filtered Molecule");
        try {
            this.filteredMoleculesPrintWriter.println("SMILES code: " + this.smilesGenerator.create(aMolecule));
        } catch (CDKException | NullPointerException anException){
            this.filteredMoleculesPrintWriter.println("SMILES code: " + this.SMILES_CODE_PLACEHOLDER);
        }
        String tmpChebiName = aMolecule.getProperty("ChEBI Name");
        if (tmpChebiName != null)
            this.filteredMoleculesPrintWriter.println("ChEBI name: " + tmpChebiName);
        String tmpChebiId = aMolecule.getProperty("ChEBI ID");
        if (tmpChebiId != null)
            this.filteredMoleculesPrintWriter.println("ChEBI ID: " + tmpChebiId);
        String tmpChemblId = aMolecule.getProperty("chembl_id");
        if (tmpChemblId != null)
            this.filteredMoleculesPrintWriter.println("ChEMBL ID: " + tmpChemblId);
        String tmpCdkTitle = aMolecule.getProperty(CDKConstants.TITLE);
        if (tmpCdkTitle != null)
            this.filteredMoleculesPrintWriter.println("CDK title: " + tmpCdkTitle);
        this.filteredMoleculesPrintWriter.println("Cause: " + aCause);
        this.filteredMoleculesPrintWriter.flush();
    }
    
    /**
     * Logs molecules that raised exceptions somewhere in the processing to the exceptions log file with exception 
     * message and stack trace, SMILES code ChEBI name and id or ChEMBL id.
     * 
     * @param anException the exception caused by the molecule
     * @param aSettingsKey a string representation of the settings tested in the current iteration, 
     * e.g. the aromaticity model
     * @param aMolecule the exception-causing molecule
     */
    private void logException(Exception anException, String aSettingsKey, IAtomContainer aMolecule) {
        this.exceptionsCounter++;
        this.exceptionsPrintWriter.println();
        this.exceptionsPrintWriter.println(this.exceptionsCounter + ". " + anException.getClass() + ": " 
                + anException.getLocalizedMessage());
        this.exceptionsPrintWriter.println("Settings key: " + aSettingsKey);
        try {
            this.exceptionsPrintWriter.println("SMILES code: " + this.smilesGenerator.create(aMolecule));
        } catch (CDKException | NullPointerException aNewException){
            this.exceptionsPrintWriter.println("SMILES code: " + this.SMILES_CODE_PLACEHOLDER);
        }
        //aMolecule.getProperty(CDKConstants.TITLE);?
        String tmpChebiName = aMolecule.getProperty("ChEBI Name");
        if (tmpChebiName != null)
            this.exceptionsPrintWriter.println("ChEBI name: " + tmpChebiName);
        String tmpChebiId = aMolecule.getProperty("ChEBI ID");
        if (tmpChebiId != null)
            this.exceptionsPrintWriter.println("ChEBI ID: " + tmpChebiId);
        String tmpChemblId = aMolecule.getProperty("chembl_id");
        if (tmpChemblId != null)
            this.exceptionsPrintWriter.println("ChEMBL ID: " + tmpChemblId);
        String tmpCdkTitle = aMolecule.getProperty(CDKConstants.TITLE);
        if (tmpCdkTitle != null)
            this.exceptionsPrintWriter.println("CDK title: " + tmpCdkTitle);
        anException.printStackTrace(this.exceptionsPrintWriter);
        this.exceptionsPrintWriter.flush();
    }
    //</editor-fold>
}

//<editor-fold defaultstate="collapsed" desc="Enum CustomAtomEncoder">
/**
 * Custom Enumeration of atom encoders for seeding atomic hash codes.
 *
 * @author Jonas Schaub
 * @see BasicAtomEncoder
 * @see AtomEncoder
 */
enum CustomAtomEncoder implements AtomEncoder {
    
    /**
     * Encode if an atom is aromatic or not. This specification is necessary to distinguish functional groups with
     * aromatic environments and those without. For example: [H]O[C] and [H]OC* (pseudo SMILES codes) should be
     * assigned different hash codes by the MoleculeHashGenerator.
     *
     * @see IAtom#isAromatic()
     */
    AROMATICITY {
        
        /**
         *{@inheritDoc}
         */
        @Override
        public int encode(IAtom anAtom, IAtomContainer aContainer) {
            return anAtom.isAromatic()? 3 : 2;
        }
    };
}
//</editor-fold>
