// =============================================
// McStasMonitorToROOT: Reads McStas monitor files into a TTree
// and provides TH2D histograms on demand.
// =============================================

#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <TSystem.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <memory>

class McStasMonitorToROOT {
public:
    // Constructor
    McStasMonitorToROOT() : fTree(nullptr), fFile(nullptr) {}

    // Destructor
    ~McStasMonitorToROOT() {
        if (fTree) delete fTree;
        if (fFile) delete fFile;
    }

    // ----------------------------
    // Add a McStas monitor file (ASCII)
    // ----------------------------
    void AddMonitorFile(const TString& filename) {
        fFilenames.push_back(filename);
    }

    // ----------------------------
    // Load all monitor files into a TTree
    // ----------------------------
    void LoadToTTree(const TString& tree_name = "monitors", const TString& file_name = "") {
        if (!fFilenames.empty() && fTree) {
            std::cerr << "Warning: TTree already loaded. Clearing existing data." << std::endl;
            Clear();
        }

        // Create a temporary ROOT file if no output file is specified
        bool use_temp_file = file_name.IsNull() || file_name == "";
        if (use_temp_file) {
            fFile = new TFile("temp_monitors.root", "RECREATE");
        } else {
            fFile = new TFile(file_name, "RECREATE");
        }

        if (!fFile->IsOpen()) {
            std::cerr << "Error: Failed to open ROOT file for TTree." << std::endl;
            return;
        }

        // Create TTree
        fTree = new TTree(tree_name, "McStas Monitor Data");

        // Define branches
        int file_id;
        double x, y, weight;
        fTree->Branch("file_id", &file_id, "file_id/I");
        fTree->Branch("x", &x, "x/D");
        fTree->Branch("y", &y, "y/D");
        fTree->Branch("weight", &weight, "weight/D");

        // Process each file
        for (size_t i = 0; i < fFilenames.size(); ++i) {
            const TString& filename = fFilenames[i];
            std::ifstream file(filename.Data());
            if (!file.is_open()) {
                std::cerr << "Error: Failed to open file: " << filename << std::endl;
                continue;
            }

            std::string line;
            while (std::getline(file, line)) {
                if (line.empty() || line[0] == '#') continue; // Skip comments/empty lines

                std::istringstream iss(line);
                if (iss >> x >> y >> weight) {
                    file_id = i; // Identify which file this entry came from
                    fTree->Fill();
                }
            }
            file.close();
        }

        // Write TTree to file
        fTree->Write();
        if (use_temp_file) {
            fFile->Close();
            // Reopen in read mode to keep the file alive
            delete fFile;
            fFile = new TFile("temp_monitors.root", "READ");
            fTree = (TTree*)fFile->Get(tree_name);
        }
    }

    // ----------------------------
    // Get a vector of TH2D histograms from the TTree
    // ----------------------------
    std::vector<TH2D*> GetTH2DHistograms(
        int n_bins_x = 100,
        double x_min = 0.0,
        double x_max = 1.0,
        int n_bins_y = 100,
        double y_min = 0.0,
        double y_max = 1.0,
        const TString& x_title = "X [m]",
        const TString& y_title = "Y [m]"
    ) {
        if (!fTree) {
            std::cerr << "Error: TTree not loaded. Call LoadToTTree first." << std::endl;
            return {};
        }

        // Clear existing histograms
        for (TH2D* hist : fHistograms) {
            delete hist;
        }
        fHistograms.clear();

        // Create a histogram for each file
        for (size_t i = 0; i < fFilenames.size(); ++i) {
            TString hist_name = Form("monitor_%zu", i);
            TString hist_title = Form("Monitor: %s", fFilenames[i].Data());

            TH2D* hist = new TH2D(
                hist_name, hist_title,
                n_bins_x, x_min, x_max,
                n_bins_y, y_min, y_max
            );
            hist->GetXaxis()->SetTitle(x_title);
            hist->GetYaxis()->SetTitle(y_title);

            fHistograms.push_back(hist);
        }

        // Set branch addresses
        int file_id;
        double x, y, weight;
        fTree->SetBranchAddress("file_id", &file_id);
        fTree->SetBranchAddress("x", &x);
        fTree->SetBranchAddress("y", &y);
        fTree->SetBranchAddress("weight", &weight);

        // Loop over all entries in the TTree
        for (Long64_t entry = 0; entry < fTree->GetEntries(); ++entry) {
            fTree->GetEntry(entry);
            if (file_id >= 0 && file_id < (int)fHistograms.size()) {
                fHistograms[file_id]->Fill(x, y, weight);
            }
        }

        return fHistograms;
    }

    // ----------------------------
    // Save histograms to a ROOT file
    // ----------------------------
    void SaveHistogramsToFile(const TString& output_filename) {
        if (fHistograms.empty()) {
            std::cerr << "Error: No histograms to save. Call GetTH2DHistograms first." << std::endl;
            return;
        }

        TFile* output_file = new TFile(output_filename, "RECREATE");
        if (!output_file->IsOpen()) {
            std::cerr << "Error: Failed to open output file: " << output_filename << std::endl;
            delete output_file;
            return;
        }

        for (TH2D* hist : fHistograms) {
            hist->Write();
        }

        output_file->Close();
        delete output_file;
        std::cout << "Histograms saved to: " << output_filename << std::endl;
    }

    // ----------------------------
    // Clear all data
    // ----------------------------
    void Clear() {
        for (TH2D* hist : fHistograms) {
            delete hist;
        }
        fHistograms.clear();

        if (fTree) delete fTree;
        fTree = nullptr;

        if (fFile) delete fFile;
        fFile = nullptr;
    }

private:
    std::vector<TString> fFilenames;  // List of monitor files
    std::vector<TH2D*> fHistograms;   // List of ROOT histograms
    TTree* fTree;                     // TTree storing monitor data
    TFile* fFile;                     // ROOT file (if used)
}

//Example of Use
void ExampleUsage() {
    // Create converter
    McStasMonitorToROOT converter;

    // Add McStas monitor files
    converter.AddMonitorFile("monitor1.dat");
    converter.AddMonitorFile("monitor2.dat");
    converter.AddMonitorFile("monitor3.dat");

    // Load data into a TTree
    converter.LoadToTTree("monitors", "monitors_tree.root");

    // Get TH2D histograms (auto-binning)
    std::vector<TH2D*> histograms = converter.GetTH2DHistograms();

    // Or specify custom binning
    // std::vector<TH2D*> histograms = converter.GetTH2DHistograms(
    //     50, -0.1, 0.1,  // X: 50 bins from -0.1 to 0.1
    //     50, -0.1, 0.1,  // Y: 50 bins from -0.1 to 0.1
    //     "X [m]", "Y [m]"
    // );

    // Save histograms to a ROOT file
    converter.SaveHistogramsToFile("monitors_histograms.root");

    // Access individual histograms
    for (TH2D* hist : histograms) {
        std::cout << "Histogram: " << hist->GetName()
                  << " (entries: " << hist->GetEntries() << ")" << std::endl;
    }
};
