/*
 * To do: Add CDK license header
 */
package org.openscience.cdk.tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
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
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.Atom;
import org.openscience.cdk.CDKTestCase;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
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
    private final String CHEBI_SD_FILE_NAME = "ChEBI_lite_3star.sdf";
    
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
     * Name of file for logging occured exceptions with causing molecules
     */
    private final String EXCEPTIONS_LOG_FILE_NAME = "Exceptions_Log";
    
    /**
     * File type of exceptions log file
     */
    private final String EXCEPTIONS_LOG_FILE_TYPE = ".txt";
    
    /**
     * First line in the exceptions log file
     */
    private final String EXCEPTIONS_LOG_FILE_HEADER = "Following molecules led to the specified exceptions:\n";
    
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
     * Seperator for file name segments (test identifier, file name, time stamp)
     */
    private final String FILE_NAME_ADDITION_SEPERATOR = "_";
    
    /**
     * Key for the master HashMap's inner maps under which to store the IAtomContainer of a functional group
     */
    private final String ATOMCONTAINER_KEY = "molecule";
    
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
     * This string will be added to an original settings key when applying the ErtlFunctionalGroupsFinder with activated 
     * generalization
     */
    private final String GENERALIZATION_SETTINGS_KEY_ADDITION = "Generalized";
    
    /**
     * Seperator for the output file's values
     */
    private final String OUTPUT_FILE_SEPERATOR = ",";
    
    /**
     * Placeholder String for every functional group's SMILES code whose real SMILES representation could not be 
     * generated
     */
    private final String SMILES_CODE_PLACEHOLDER = "SMILES code could not be created.";
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
     * Pseudo SMILES representation of an undifined pseudo atom
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
     * True if the filtered molecules were logged in the filtered molecules log file; This is only neccessary in the 
     * first iteration since the applied filters are the same in every iteration (assuming that in a single test
     * only one SD file is analyzed)
     */
    private boolean areFilteredMoleculesLogged;
    
    /**
     * True if all operations in initialize() were successfull and the test is able to run
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
     * Instance of the ErtlFunctionalGroupsFinder
     */
    private ErtlFunctionalGroupsFinder ertlFGFinder;
    
    /**
     * Master HashMap for storing results; Its keys are the hash codes produced by the MoleculeHashGenerator for the 
     * functional groups and its values are inner HashMaps that hold the IAtomContainer of a functional group and its 
     * frequencies for different settings as String-Object pairs
     */
    private HashMap<Long, HashMap> masterHashMap;
    
    /**
     * A list for storing all used settings keys in a test
     */
    private Set<String> settingsKeysList;
    
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
    /**
     * Test for analyzing ChEBI database for all four different aromaticity models supplied by the cdk: daylight, cdk, 
     * piBonds and cdkAllowingExocyclic
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     * 
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occures
     */
    @Ignore
    @Test
    public void analyzeChebi() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CHEBI_TEST_IDENTIFIER);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file named: " + this.CHEBI_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEBI_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 4;
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
        
        boolean tmpAreMultiplesCounted = true;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Test for analyzing ChEBI database for all four different aromaticity models supplied by the cdk: daylight, cdk, 
     * piBonds and cdkAllowingExocyclic;
     * <p>
     * Difference to analyzeChebi(): If the same functional group occurres multiple times in the same molecule 
     * it is counted only once
     * 
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occures
     */
    @Test
    public void analyzeChebiNoMultiples() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CHEBI_NO_MULTIPLES_TEST_IDENTIFIER);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file named: " + this.CHEBI_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEBI_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 4;
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
        
        boolean tmpAreMultiplesCounted = false;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        
        System.out.println("\nAll analyses are done!");
        this.postProcessAndSaveData();
        System.out.println("\nFinished!");
        System.out.println("\nNumber of occured exceptions: " + this.exceptionsCounter);
    }
    
    /**
     * Test for analyzing ChEMBL database for all four different aromaticity models supplied by the cdk: daylight, cdk, 
     * piBonds and cdkAllowingExocyclic
     * <p>
     * (Functional groups occurring multiple times in the same molecule are counted multiple times)
     * 
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occures
     */
    @Ignore
    @Test
    public void analyzeChembl() throws Exception {
        this.initialize(this.CHEMBL_SD_FILE_NAME, this.CHEMBL_TEST_IDENTIFIER);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file named: " + this.CHEMBL_SD_FILE_NAME);
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpChebiSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + this.CHEMBL_SD_FILE_NAME)
                .getFile());
        int tmpRequiredNumberOfReaders = 4;
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
        
        boolean tmpAreMultiplesCounted = true;
        
        this.calculateAbsoluteFGFrequencies(tmpReaders[0], this.DAYLIGHT_MODEL_SETTINGS_KEY, tmpDaylightModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[1], this.CDK_MODEL_SETTINGS_KEY, tmpCdkModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[2], this.PIBONDS_MODEL_SETTINGS_KEY, tmpPiBondsModel, tmpAreMultiplesCounted);
        this.calculateAbsoluteFGFrequencies(tmpReaders[3], this.CDK_EXOCYCLIC_MODEL_SETTINGS_KEY, tmpCdkAllowingExocyclicModel, tmpAreMultiplesCounted);
        
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
     * @throws java.lang.Exception if initialize() throws an exception or an unexpected exception occures
     */
    @Ignore
    @Test
    public void cycleFinderTest() throws Exception {
        this.initialize(this.CHEBI_SD_FILE_NAME, this.CYCLE_FINDER_TEST_IDENTIFIER);
        Assume.assumeTrue(this.isTestAbleToRun);
        
        System.out.println("\nLoading file named: " + this.CHEBI_SD_FILE_NAME);
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
    
    //<editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Initializes all class variables and determines the output directory
     * 
     * @param aFileName name of the SD file to analyze for a quick pre-check if it is present and the test is therefore 
     * meant to run
     * @param aTestIdentifier a folder with this name will be created in the output directory and it will be added to 
     * the output and log files' names for association of test and files
     * @throws java.lang.Exception if one the FileWriter instances can not be instantiated or more than 
     * Integer.MAX-VALUE tests are to be run this minute (error in the naming of output files) or an unexpected 
     * exception occures.
     */
    private void initialize(String aFileName, String aTestIdentifier) throws Exception {
        System.out.println("\n#########################################################################\n");
        System.out.println("Starting new test, identifier: " + aTestIdentifier);
        System.out.println("\nInitializing class variables...");
        this.isTestAbleToRun = true;
        //First, check if the SD file is present and ignore test if it is not
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = new File(tmpClassLoader.getResource(this.SDF_RESOURCES_PATH + aFileName).getFile());
        if (!tmpSDFile.exists()) {
            System.out.println("\n\tUnable to find a file named " + aFileName 
                    + " in the resources directory. Test is ignored.");
            this.isTestAbleToRun = false;
            Assume.assumeTrue(false);
            return;
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
        System.out.println("\n\tExceptions will be written to: " + tmpExceptionsLogFile.getName());
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
        //examplary fast setup of the MoleculeHashGenerator
        this.molHashGenerator = new HashGeneratorMaker()
                .depth(8)
                .elemental()
                .isotopic()
                .charged()
                .orbital()
                .molecular();
        //Generalization is later done by invoking ErtlFunctionalGroupsFinder.expandGeneralizedEnvironments()
        this.ertlFGFinder = new ErtlFunctionalGroupsFinder(Mode.NO_GENERALIZATION);
        this.masterHashMap = new HashMap(this.MASTER_HASHMAP_INITIAL_CAPACITY, this.MASTER_HASHMAP_LOAD_FACTOR);
        this.settingsKeysList = new HashSet<>();
        this.exceptionsCounter = 0;
        this.areFilteredMoleculesLogged = false;
        String[] tmpMetalNumbersStrings = this.NON_METALLIC_ATOMIC_NUMBERS.split(",");
        Integer[] tmpMetalNumbersInt = new Integer[tmpMetalNumbersStrings.length];
        for (int i = 0; i < tmpMetalNumbersStrings.length; i++) {
            tmpMetalNumbersInt[i] = Integer.parseInt(tmpMetalNumbersStrings[i]);
        }
        this.nonMetallicAtomicNumbersSet = new HashSet(Arrays.asList(tmpMetalNumbersInt));
        System.out.println("\n\tDone initializing class variables.");
    }
    
    /**
     * Does one iteration of loading molecules from an IteratingSDFReader, applying the given aromaticity model, 
     * extracting their functional groups (with and without generalization) and adding the results to the master HashMap;
     * Exceptions caused by read-in molecules are caught and logged
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
        int tmpMoleculesCounter = 0; //total number of analyzed molecules
        int tmpChargeCounter = 0; //number of molecules that contain one or more charged atoms
        int tmpUnconnectedCounter = 0; //number of structures with unconnected substructures
        int tmpFilteredMoleculesCounter = 0; //filtered molecules, see below
        int tmpMetallicCounter = 0; //number of molecules that were filtered because they contain one or more metal, metalloid or "R" atoms
        int tmpNoAtomOrBondCounter = 0; //number of molecules that were filtered because they do not contain an atom or a bond
        int tmpSkippedMoleculesCounter = 0; //molecules that caused exceptions and are therefore skipped
        int tmpNoFunctionalGroupsCounter = 0; //number of molecules that contain no functional group
        //Note: A molecule can be supposed to be filtered for more than one of the named reasons but only the first 
        //tested reason will be named as cause
        while (aReader.hasNext()) {
            tmpMoleculesCounter++;
            List<IAtomContainer> tmpFunctionalGroups;
            IAtomContainer tmpMolecule = (IAtomContainer) aReader.next();
            if (tmpMolecule == null) {
                break;
            }
            String tmpCause = "";
            IAtomContainer tmpOriginalMolecule = null;
            boolean tmpIsBiggestFragmentSelected = false;
            try {
                //Filter molecules with atom or bond count zero
                if (tmpMolecule.getAtomCount() == 0 || tmpMolecule.getBondCount() == 0) {
                    tmpCause = "Atom or bond count 0";
                    tmpFilteredMoleculesCounter++;
                    tmpNoAtomOrBondCounter++;
                    if (!this.areFilteredMoleculesLogged) {
                        this.logFilteredMolecule(tmpMolecule, tmpFilteredMoleculesCounter, tmpCause);
                    }
                    continue;
                }
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMolecule);
                //Preprocessing: From structures containing two or more unconnected structures (e.g. ions) 
                //choose the largest structure for analysis
                //This step must be done prior to the next step!
                if (!ConnectivityChecker.isConnected(tmpMolecule)) {
                    tmpUnconnectedCounter++;
                    IAtomContainerSet tmpFragmentsSet = ConnectivityChecker.partitionIntoMolecules(tmpMolecule);
                    IAtomContainer tmpBiggestFragment = null;
                    for (IAtomContainer tmpFragment : tmpFragmentsSet.atomContainers()) {
                        if (tmpBiggestFragment == null || tmpBiggestFragment.getAtomCount() < tmpFragment.getAtomCount()) {
                            tmpBiggestFragment = tmpFragment;
                        }
                    }
                    tmpIsBiggestFragmentSelected = true;
                    tmpOriginalMolecule = tmpMolecule.clone();
                    tmpMolecule = tmpBiggestFragment;
                }
                //Filter molecules containing metals, metalloids or "R" atoms
                boolean tmpIsForbiddenAtomicNumberDetected = false;
                for (IAtom tmpAtom : tmpMolecule.atoms()) {
                    if (!this.nonMetallicAtomicNumbersSet.contains(tmpAtom.getAtomicNumber())) {
                        tmpIsForbiddenAtomicNumberDetected = true;
                        break;
                    }
                }
                if (tmpIsForbiddenAtomicNumberDetected) {
                    tmpCause = "Contains one or more metal, metalloid or \"R\" atoms";
                    tmpFilteredMoleculesCounter++;
                    tmpMetallicCounter++;
                    if (!this.areFilteredMoleculesLogged) {
                        if (tmpIsBiggestFragmentSelected && tmpOriginalMolecule != null) {
                            this.logFilteredMolecule(tmpOriginalMolecule, tmpFilteredMoleculesCounter, tmpCause);
                        } else {
                            this.logFilteredMolecule(tmpMolecule, tmpFilteredMoleculesCounter, tmpCause);
                        }
                    }
                    continue;
                }
                //Preprocessing: Neutralize charges
                for (IAtom tmpAtom : tmpMolecule.atoms()) {
                    if (tmpAtom.getFormalCharge() != 0) {
                        tmpChargeCounter++;
                        tmpAtom.setFormalCharge(0);
                        CDKHydrogenAdder tmpHAdder = CDKHydrogenAdder.getInstance(tmpMolecule.getBuilder());
                        CDKAtomTypeMatcher tmpMatcher = CDKAtomTypeMatcher.getInstance(tmpMolecule.getBuilder());
                        IAtomType tmpMatchedType = tmpMatcher.findMatchingAtomType(tmpMolecule, tmpAtom);
                        if (tmpMatchedType != null) {
                            AtomTypeManipulator.configure(tmpAtom, tmpMatchedType);
                        }
                        tmpHAdder.addImplicitHydrogens(tmpMolecule, tmpAtom);
                    }
                }
                //Now the analysis of functional groups:
                anAromaticity.apply(tmpMolecule);
                this.ertlFGFinder.find(tmpMolecule);
                this.ertlFGFinder.markAtoms(tmpMolecule);
                tmpFunctionalGroups = this.ertlFGFinder.extractGroups(tmpMolecule);
            } catch (Exception anException) {
                tmpSkippedMoleculesCounter++;
                if (tmpIsBiggestFragmentSelected && tmpOriginalMolecule != null) {
                    this.logException(anException, aSettingsKey, tmpOriginalMolecule);
                } else {
                    this.logException(anException, aSettingsKey, tmpMolecule);
                }
                continue;
            }
            if (tmpFunctionalGroups == null || tmpFunctionalGroups.isEmpty()) {
                tmpNoFunctionalGroupsCounter++;
                continue;
            }
            this.addFunctionalGroupsToMasterMap(tmpFunctionalGroups, aSettingsKey, anAreMultiplesCounted);
            //Generalize environments of the extracted functional groups
            List<IAtomContainer> tmpFunctionalGroupsGeneralized = new LinkedList<>();
            for (IAtomContainer tmpFunctionalGroup : tmpFunctionalGroups) {
                try {
                    tmpFunctionalGroupsGeneralized.add(tmpFunctionalGroup.clone());
                } catch (CloneNotSupportedException aCloneNotSupportedException){
                    if (tmpIsBiggestFragmentSelected && tmpOriginalMolecule != null) {
                        this.logException(aCloneNotSupportedException, 
                                aSettingsKey + " (Cloning one of the molecule's functional groups)", 
                                tmpOriginalMolecule);
                    } else {
                        this.logException(aCloneNotSupportedException, 
                                aSettingsKey + " (Cloning one of the molecule's functional groups)", 
                                tmpMolecule);
                    }
                }
            }
            try {
                tmpFunctionalGroupsGeneralized = this.ertlFGFinder.expandGeneralizedEnvironments(tmpFunctionalGroupsGeneralized);
            } catch (NullPointerException anException) {
                if (tmpIsBiggestFragmentSelected && tmpOriginalMolecule != null) {
                    this.logException(anException, 
                            aSettingsKey + " (Generalizing the molecule's functional groups)", 
                            tmpOriginalMolecule);
                } else {
                    this.logException(anException, 
                            aSettingsKey + " (Generalizing the molecule's functional groups)", 
                            tmpMolecule);
                }
            }
            this.addFunctionalGroupsToMasterMap(tmpFunctionalGroupsGeneralized, 
                    aSettingsKey + this.GENERALIZATION_SETTINGS_KEY_ADDITION,
                    anAreMultiplesCounted);
        }
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
     * Inserts a list of IAtomContainers (the functional groups of one molecule) into the master HashMap; If the 
     * functional group is already inserted (with this settings key) its frequency for the given settings key is raised 
     * by one or else a new inner HashMap will be created for it
     * 
     * @param aFunctionalGroupsList the functional groups of one molecule to be inserted
     * @param aSettingsKey will be the key of the inner HashMap inside the master HashMap
     * @param anAreMultiplesCounted if false, functional groups that occur multiple times in aFunctionalGroupsList will 
     * only be entered once into the master Hashmap
     */
    private void addFunctionalGroupsToMasterMap(
            List<IAtomContainer> aFunctionalGroupsList, 
            String aSettingsKey, 
            boolean anAreMultiplesCounted) {
        
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
                HashMap tmpNewInnerMap = new HashMap(this.INNER_HASHMAPS_INITIAL_CAPACITY);
                tmpNewInnerMap.put(this.ATOMCONTAINER_KEY, tmpFunctionalGroup);
                tmpNewInnerMap.put(aSettingsKey, 1);
                this.masterHashMap.put(tmpHashCode, tmpNewInnerMap);
            }
            tmpAlreadyEnteredFGsForThisMol.add(tmpHashCode);
        }
        
    }
    
    /**
     * Writes all frequency data with the respective hash code, SMILES code and pseudo SMILES code for all functional 
     * groups in the master HashMap to the output file
     * <p>
     * Note: The IAtomContainer object stored in the master HashMap's inner maps are altered in this method for pseudo 
     * SMILES creation! And all PrintWriter instances will be closed
     */
    private void postProcessAndSaveData() {
        System.out.println("\nWriting to file...");
        //Writing the output file's header
        String tmpFileHeader = this.HASH_CODE_KEY + this.OUTPUT_FILE_SEPERATOR + this.PSEUDO_SMILES_CODE_KEY 
                + this.OUTPUT_FILE_SEPERATOR + this.SMILES_CODE_KEY;
        for (String tmpSettingsKey : this.settingsKeysList) {
            tmpFileHeader += this.OUTPUT_FILE_SEPERATOR + tmpSettingsKey;
        }
        this.dataOutputPrintWriter.println(tmpFileHeader);
        this.dataOutputPrintWriter.flush();
        Iterator tmpFunctionalGroupsIterator = this.masterHashMap.keySet().iterator();
        //Iteration for all molecules in the master HashMap
        while (tmpFunctionalGroupsIterator.hasNext()) {
            long tmpHashCode = (long)tmpFunctionalGroupsIterator.next();
            HashMap tmpInnerMap = (HashMap)this.masterHashMap.get(tmpHashCode);
            String tmpSmilesCode;
            String tmpPseudoSmilesCode;
            IAtomContainer tmpFunctionalGroup = (IAtomContainer)tmpInnerMap.get(this.ATOMCONTAINER_KEY);
            try {
                //Creation of unique SMILES code
                tmpSmilesCode = this.smilesGenerator.create(tmpFunctionalGroup);
                //Creation of pseudo SMILES code
                for (IAtom tmpAtom : tmpFunctionalGroup.atoms()) {
                    if (tmpAtom.isAromatic()) {
                        IAtom tmpReplacementAtom = null;
                        if ("C".equals(tmpAtom.getSymbol())) {
                            tmpReplacementAtom = new Atom("Ce");
                        } else if ("N".equals(tmpAtom.getSymbol())) {
                            tmpReplacementAtom = new Atom("Nd");
                        } else if ("S".equals(tmpAtom.getSymbol())) {
                            tmpReplacementAtom = new Atom("Sm");
                        } else if ("O".equals(tmpAtom.getSymbol())) {
                            tmpReplacementAtom = new Atom("Os");
                        }
                        if (tmpReplacementAtom != null) {
                            Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                            AtomContainerManipulator.replaceAtomByAtom(tmpFunctionalGroup, tmpAtom, tmpReplacementAtom);
                            tmpReplacementAtom.setImplicitHydrogenCount(
                                    tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                        }
                    }
                    tmpAtom.setIsAromatic(false);
                    if (tmpAtom instanceof IPseudoAtom) {
                        IAtom tmpReplacementAtom = new Atom("Es");
                        Integer tmpImplicitHydrogenCount = tmpAtom.getImplicitHydrogenCount();
                        AtomContainerManipulator.replaceAtomByAtom(tmpFunctionalGroup, tmpAtom, tmpReplacementAtom);
                        tmpReplacementAtom.setImplicitHydrogenCount(
                                tmpImplicitHydrogenCount == null ? 0 : tmpImplicitHydrogenCount);
                    }
                }
                for (IBond tmpBond : tmpFunctionalGroup.bonds()) {
                    tmpBond.setIsAromatic(false);
                }
                tmpPseudoSmilesCode = this.smilesGenerator.create(tmpFunctionalGroup);
                tmpPseudoSmilesCode = tmpPseudoSmilesCode
                        .replaceAll("(\\[Es\\])", this.PSEUDO_SMILES_R_ATOM)
                        .replaceAll("(\\[Ce\\])", this.PSEUDO_SMILES_AROMATIC_CARBON)
                        .replaceAll("(\\[Nd\\])", this.PSEUDO_SMILES_AROMATIC_NITROGEN)
                        .replaceAll("(\\[Sm\\])", this.PSEUDO_SMILES_AROMATIC_SULPHUR)
                        .replaceAll("(\\[Os\\])", this.PSEUDO_SMILES_AROMATIC_OXYGEN);
            } catch (CDKException | NullPointerException anException) {
                this.logException(anException, "Creating SMILES code, no settings", tmpFunctionalGroup);
                tmpSmilesCode = this.SMILES_CODE_PLACEHOLDER;
                tmpPseudoSmilesCode = this.SMILES_CODE_PLACEHOLDER;
            }
            //Writing the record for this functional group
            String tmpRecord = tmpHashCode + this.OUTPUT_FILE_SEPERATOR + tmpPseudoSmilesCode 
                    + this.OUTPUT_FILE_SEPERATOR + tmpSmilesCode;
            for (String tmpSettingsKey : this.settingsKeysList) {
                if (tmpInnerMap.get(tmpSettingsKey) == null) {
                    tmpInnerMap.put(tmpSettingsKey, 0);
                }
                tmpRecord += this.OUTPUT_FILE_SEPERATOR + tmpInnerMap.get(tmpSettingsKey);
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
     * Logs molecules that are filtered from the SD file to the filtered molecules file with SMILES code, ChEBI name 
     * and id or ChEMBL id and why they were filtered
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
            this.filteredMoleculesPrintWriter.println("SMILES code: " +this.SMILES_CODE_PLACEHOLDER);
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
        this.filteredMoleculesPrintWriter.println("Cause: " + aCause);
        this.filteredMoleculesPrintWriter.flush();
    }
    
    /**
     * Logs molecules that raised exceptions somewhere in the processing to the exceptions log file with exception 
     * message and stack trace, SMILES code ChEBI name and id or ChEMBL id
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
            this.exceptionsPrintWriter.println("SMILES code: " +this.SMILES_CODE_PLACEHOLDER);
        }
        String tmpChebiName = aMolecule.getProperty("ChEBI Name");
        if (tmpChebiName != null)
            this.exceptionsPrintWriter.println("ChEBI name: " + tmpChebiName);
        String tmpChebiId = aMolecule.getProperty("ChEBI ID");
        if (tmpChebiId != null)
            this.exceptionsPrintWriter.println("ChEBI ID: " + tmpChebiId);
        String tmpChemblId = aMolecule.getProperty("chembl_id");
        if (tmpChemblId != null)
            this.exceptionsPrintWriter.println("ChEMBL ID: " + tmpChemblId);
        anException.printStackTrace(this.exceptionsPrintWriter);
        this.exceptionsPrintWriter.flush();
    }
    //</editor-fold>
}
