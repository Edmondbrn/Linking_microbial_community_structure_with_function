import java.io.BufferedInputStream;
import java.io.File;
import java.util.Scanner;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

public class Parsing {

    private static void ExtractLines(ArrayList<String> liste_files, String path) {
        // Function to extract lines from the file and then store them in a list
        int nb_files = liste_files.size(); // variable to print the progress of the parsing
        int cpt = 1;
        for (String files : liste_files){ // loop to parse each file
            HashMap<String, String> dict = new HashMap<>(); // dictionnary to do the link between various element of the file
            int i = 0;
            dict.put("VERSION", files); // name of the file or the cluster
            try{
                File file = new File(path + "/" + files); // open and read the genbank file
                Scanner scanner = new Scanner(file);
                while (scanner.hasNextLine()){
                    String line = scanner.nextLine().strip(); // remove garbage characaters at the end and at the beginning of the line
                     if (line.startsWith("ORGANISM")) // if we find the organism name, we store it in the dictionnary
                        dict.put("ORGANISM", line.substring(line.indexOf(" ") +1));   // Store the organisme name by cutting the line with the fisrt space as reference
                    
                    else if (line.startsWith("aSDomain")){ // if we find the AMP-binding domain, we store the sequence in the dictionnary
                        line = scanner.nextLine().strip();
                        if (line.startsWith("/aSDomain=\"AMP-binding\"")){
                            while(!line.startsWith("/translation")){ // we search the sequence of the AMP-binding domain (we pass lines until we find the sequence)
                                line = scanner.nextLine().strip();
                            }
                            String sequence = new String(line.split("\"")[1] + "\n"); // we only take the part after the " character
                            while(!line.endsWith("\"")){ // we take the sequence line by line until we find the end of the sequence
                                line = scanner.nextLine().strip();
                                sequence += line + "\n";
                            }
                            dict.put("AMP-Domain_" + i, sequence); // give a number to the AMP-binding domain to avoid duplicates
                            i++;
                        }
                    }
                }
                scanner.close();
            } catch (Exception e) {
                System.out.println("Error, cannot find the speciefied file :" + path + " Please check the path and try again");
                e.printStackTrace();
            }
            writeFile(dict); // write the information in the output file
            verbose(cpt, nb_files); // display the progress of the parsing
            cpt ++;
        }
    }

    private static void verbose(int cpt, int nb_files){
        // Function to display the progress of the parsing
        System.out.println("File number " + cpt + " out of " + nb_files + " has been parsed");
        if (cpt == nb_files)
            System.out.println("Parsing completed successfully");
    }

    private static void writeFile(HashMap<String, String> data_dict){
        /*
         * Function to write information in the output file
         */
        int i =0;
        try{
            File out_fil = new File("parsed_AMP_Mibig.txt");
            FileWriter writer = new FileWriter(out_fil, true);
            for (String key: data_dict.keySet()){
                if (key.startsWith("AMP-Domain")){
                    writer.write(">" + data_dict.get("VERSION") + "_cluster_" + i  + data_dict.get("ORGANISM") + "\n");
                    writer.write(data_dict.get(key) + "\n");
                    i ++;
                }
            } writer.close();
        } catch (IOException e) {
            System.out.println("An error occured while writing the file");
            e.printStackTrace();
        }
    }

    private static ArrayList<String> Extractfiles(String path){
        // Function to extract files from the directory and then store them in a list
        ArrayList<String> liste_files = new ArrayList<String>();
        try{
            File file = new File(path); // creation of an object to manipulate files
            File[] files = file.listFiles(); // list of files in the directory
            for (File f : files){
                if (f.isFile()){
                    liste_files.add(f.getName());
                }
            }
        } catch (Exception e) {
            System.out.println("Error, cannot find the speciefied directory :" + path + " Please check the path and try again");
            e.printStackTrace();
        }
        return liste_files;
}

    private static void ExtractName(){
        // Function to extract the name of the organism and to store them in a file
        String path = "parsed_AMP_antismshdb.txt";
        File file = new File(path);
        ArrayList<String> organism_list = new ArrayList<>();
        try{
            Scanner scanner = new Scanner(file);
            FileWriter writer = new FileWriter("Organism_list.txt");
            while (scanner.hasNextLine()){
                String line = scanner.nextLine().strip();
                if (line.startsWith(">")){ // extract name from the header
                    String organism = line.substring(line.indexOf(" ") +1);
                    if (!organism_list.contains(organism)){ // verify if the organism has already been added to the list
                        organism_list.add(organism);
                        writer.write(organism + "\n");
                        System.out.println(organism);
                    }
                }
            }
            scanner.close();
            writer.close();
        } catch (Exception e) {
            System.out.println("Error, cannot find the speciefied file :" + path + " Please check the path and try again");
            e.printStackTrace();
        }
        System.out.println("Organism list has been extracted successfully");
    }

    private static void LinkeOrgaTax(String path_orga, String path_id){
        /*
         * Function to create the file that will do the link between organism and taxid
         */
        File file_orga = new File(path_orga);
        File file_id = new File(path_id);
        ArrayList<String> cluster_list = new ArrayList<>();
        try{
            Scanner scanner_orga = new Scanner(file_orga);
            FileWriter writer = new FileWriter("antismashdb_taxa_dicto.txt", true);
            int i = 1;
            while (scanner_orga.hasNextLine()){
                String line_orga = scanner_orga.nextLine().strip();
                if (line_orga.startsWith(">")){
                    String organism = line_orga.substring(line_orga.indexOf(" ") +1).strip(); // get the organism name from the header
                    Scanner scanner_id = new Scanner(file_id);
                    String cluster_name = line_orga.substring(1, line_orga.indexOf(".gbk") +4).strip(); // get the cluster name from the header
                    if (cluster_list.contains(cluster_name)) // verify if the cluster has already been added to the list (to avoid duplicates)
                        continue;
                    while (scanner_id.hasNextLine()){
                        String line_id = scanner_id.nextLine().strip();
                        if (line_id.contains(organism)){
                            line_id = line_id.substring(line_id.indexOf("|")+1).strip(); // remove the three columns from NCBI ouput
                            line_id = line_id.substring(line_id.indexOf("|")+1).strip();
                            line_id = line_id.substring(line_id.indexOf("|")+1).strip();
                            writer.write(cluster_name + "\t" + "b\'" + line_id + "\'\\n\n");
                            break;
                        }
                    }
                    scanner_id.close();
                    cluster_list.add(cluster_name);
                    System.out.println("Cluster number " + i +  " has been linked to its taxid");
                    i++;
                }
            }
            scanner_orga.close();
            writer.close();
        } catch(Exception e){
            System.out.println("Error, cannot find the speciefied file :" + path_orga + " Please check the path and try again");
            e.printStackTrace();
        }

    }

    private static void productFile(String path_id){
        /*
         * Function to assign product to the cluster (today only NRPS and PKS are considered)
         */
        File file_id = new File(path_id);
        try{
            Scanner scanner = new Scanner(file_id);
            FileWriter writer = new FileWriter("asdb_parsed_product_list.txt", true);
            while (scanner.hasNextLine()){
               String line = scanner.nextLine().strip().split("\t")[0];
               writer.write(line + "\t" + "nrps" + "\n");
            }
            scanner.close();
            writer.close();
        } catch (Exception e){
            System.out.println("Error, cannot find the speciefied file :" + path_id + " Please check the path and try again");
            e.printStackTrace();
        }
        System.out.println("Product list has been created successfully");
    }


    public static void main(String[] args) {
        String path = new String("/home/edmond/database_update/Parsing_database/ressources/data/mibig_db/mibig_gbk_3.1");
        ArrayList<String> file_list = Extractfiles(path);
        ExtractLines(file_list, path);
        ExtractName();
        //LinkeOrgaTax("/home/edmond/database_update/Parsing_database/src/parsed_AMP_antismshdb.txt", "/home/edmond/database_update/Parsing_database/src/tax_report (1).txt");
        //productFile("/home/edmond/database_update/Parsing_database/src/antismashdb_taxa_dicto.txt");
    }

}
