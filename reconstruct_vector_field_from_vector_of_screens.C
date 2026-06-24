// =============================================
// SIRT Reconstruction Skeleton (ROOT C++)
// Voxels contain 3D vectors (e.g., for direction, polarization, etc.)
// =============================================

#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <memory>

// ----------------------------
// Voxel struct: Stores a 3D vector (x, y, z components)
// ----------------------------
struct Voxel {
    double value;      // Scalar value (e.g., intensity, density)
    double vector[3]; // 3D vector (e.g., direction, polarization)

    Voxel() : value(0.0) {
        vector[0] = vector[1] = vector[2] = 0.0;
    }
};

// ----------------------------
// 3D Volume class: Stores voxels with 3D vectors
// ----------------------------
class VectorVolume3D {
public:
    VectorVolume3D(int nx, double x_min, double x_max,
                    int ny, double y_min, double y_max,
                    int nz, double z_min, double z_max)
        : fNx(nx), fNy(ny), fNz(nz),
          fXMin(x_min), fXMax(x_max),
          fYMin(y_min), fYMax(y_max),
          fZMin(z_min), fZMax(z_max) {
        fVoxels.resize(nx * ny * nz);
    }

    // Get voxel at (x, y, z)
    Voxel& GetVoxel(int x, int y, int z) {
        return fVoxels[x + y * fNx + z * fNx * fNy];
    }

    // Get voxel at (x, y, z) (const version)
    const Voxel& GetVoxel(int x, int y, int z) const {
        return fVoxels[x + y * fNx + z * fNx * fNy];
    }

    // Get dimensions
    int GetNx() const { return fNx; }
    int GetNy() const { return fNy; }
    int GetNz() const { return fNz; }
    double GetXMin() const { return fXMin; }
    double GetXMax() const { return fXMax; }
    double GetYMin() const { return fYMin; }
    double GetYMax() const { return fYMax; }
    double GetZMin() const { return fZMin; }
    double GetZMax() const { return fZMax; }

    // Convert to TH3D (scalar part only)
    TH3D* ToTH3D(const TString& name, const TString& title) const {
        TH3D* hist = new TH3D(
            name, title,
            fNx, fXMin, fXMax,
            fNy, fYMin, fYMax,
            fNz, fZMin, fZMax
        );
        for (int x = 0; x < fNx; ++x) {
            for (int y = 0; y < fNy; ++y) {
                for (int z = 0; z < fNz; ++z) {
                    hist->SetBinContent(x+1, y+1, z+1, GetVoxel(x, y, z).value);
                }
            }
        }
        return hist;
    }

private:
    int fNx, fNy, fNz;
    double fXMin, fXMax, fYMin, fYMax, fZMin, fZMax;
    std::vector<Voxel> fVoxels;
};

// ----------------------------
// SIRT Reconstruction Class
// ----------------------------
class SIRTReconstruction {
public:
    SIRTReconstruction(int nx, double x_min, double x_max,
                        int ny, double y_min, double y_max,
                        int nz, double z_min, double z_max)
        : fVolume(nx, x_min, x_max, ny, y_min, y_max, nz, z_min, z_max) {}

    // ----------------------------
    // Set the initial volume (3D vector field)
    // ----------------------------
    void SetInitialVolume(const VectorVolume3D& volume) {
        fVolume = volume;
    }

    // ----------------------------
    // Custom Forward Projection:
    // Maps 3D voxels (with 3D vectors) to 2D pixels on multiple screens.
    // Replace this with your actual forward projection logic.
    // ----------------------------
    std::vector<TH2D*> ForwardProject(const std::vector<std::pair<int, int>>& screen_dimensions) {
        std::vector<TH2D*> projections;
        for (size_t screen_idx = 0; screen_idx < screen_dimensions.size(); ++screen_idx) {
            int n_pixels_x = screen_dimensions[screen_idx].first;
            int n_pixels_y = screen_dimensions[screen_idx].second;

            TH2D* projection = new TH2D(
                Form("projection_screen_%zu", screen_idx),
                Form("Forward Projection (Screen %zu)", screen_idx),
                n_pixels_x, 0, n_pixels_x,
                n_pixels_y, 0, n_pixels_y
            );
            projections.push_back(projection);

            // Loop over all voxels
            for (int x = 0; x < fVolume.GetNx(); ++x) {
                for (int y = 0; y < fVolume.GetNy(); ++y) {
                    for (int z = 0; z < fVolume.GetNz(); ++z) {
                        const Voxel& voxel = fVolume.GetVoxel(x, y, z);

                        // TODO: Replace with your custom voxel-to-pixel mapping logic.
                        // Example: Use the 3D vector in the voxel and screen geometry to determine pixel contributions.
                        // Here, we use a simplified projection for demonstration.

                        // Example: Project voxel center to pixel coordinates for this screen.
                        double pixel_x = (x + 0.5) * (n_pixels_x / fVolume.GetNx());
                        double pixel_y = (y + 0.5) * (n_pixels_y / fVolume.GetNy());

                        // Clamp to pixel bounds
                        int px = TMath::Min(TMath::Max(1, (int)(pixel_x + 0.5)), n_pixels_x);
                        int py = TMath::Min(TMath::Max(1, (int)(pixel_y + 0.5)), n_pixels_y);

                        // Accumulate voxel contribution to the pixel (scalar part)
                        projection->AddBinContent(px, py, voxel.value);

                        // TODO: Optionally, use the 3D vector to modify the contribution.
                        // Example: projection->AddBinContent(px, py, voxel.value * voxel.vector[0]);
                    }
                }
            }
        }
        return projections;
    }

    // ----------------------------
    // Backward Projection (Adjoint):
    // Maps errors from all screens back to 3D voxels.
    // Replace this with the adjoint of your forward projection.
    // ----------------------------
    void BackwardProject(const std::vector<TH2D*>& errors, VectorVolume3D& volume) {
        // Loop over all screens
        for (size_t screen_idx = 0; screen_idx < errors.size(); ++screen_idx) {
            const TH2D* error = errors[screen_idx];

            // Loop over all pixels in this screen
            for (int px = 0; px < error->GetNbinsX(); ++px) {
                for (int py = 0; py < error->GetNbinsY(); ++py) {
                    double pixel_error = error->GetBinContent(px+1, py+1);

                    // TODO: Replace with your custom pixel-to-voxel mapping logic.
                    // Example: Distribute the pixel error to all voxels that contributed to it.
                    // Here, we use a simplified backprojection for demonstration.

                    // Example: Map pixel back to voxels (uniformly)
                    for (int x = 0; x < volume.GetNx(); ++x) {
                        for (int y = 0; y < volume.GetNy(); ++y) {
                            for (int z = 0; z < volume.GetNz(); ++z) {
                                // TODO: Replace with your actual adjoint logic.
                                // Example: Weight the contribution based on your forward model.
                                double weight = 1.0 / (volume.GetNx() * volume.GetNy() * volume.GetNz());
                                volume.GetVoxel(x, y, z).value += pixel_error * weight;

                                // TODO: Optionally, update the 3D vector components.
                                // Example: volume.GetVoxel(x, y, z).vector[0] += pixel_error * weight * some_factor;
                            }
                        }
                    }
                }
            }
        }
    }

    // ----------------------------
    // SIRT Reconstruction
    // ----------------------------
    VectorVolume3D SIRTReconstruct(
        const std::vector<TH2D*>& measured_data, // Measured pixel data (one TH2D per screen)
        int max_iterations = 10,
        double relaxation = 0.1                  // Relaxation parameter (0 < λ ≤ 1)
    ) {
        // Clone initial volume
        VectorVolume3D reconstructed = fVolume;

        // Precompute normalization (sum of forward projection weights for each voxel across all screens)
        VectorVolume3D normalization(fVolume.GetNx(), fVolume.GetXMin(), fVolume.GetXMax(),
                                      fVolume.GetNy(), fVolume.GetYMin(), fVolume.GetYMax(),
                                      fVolume.GetNz(), fVolume.GetZMin(), fVolume.GetZMax());

        // Create a "ones" volume for normalization
        VectorVolume3D ones(fVolume.GetNx(), fVolume.GetXMin(), fVolume.GetXMax(),
                            fVolume.GetNy(), fVolume.GetYMin(), fVolume.GetYMax(),
                            fVolume.GetNz(), fVolume.GetZMin(), fVolume.GetZMax());
        for (int x = 0; x < ones.GetNx(); ++x) {
            for (int y = 0; y < ones.GetNy(); ++y) {
                for (int z = 0; z < ones.GetNz(); ++z) {
                    ones.GetVoxel(x, y, z).value = 1.0;
                }
            }
        }

        // Forward project "ones" to all screens to get normalization weights
        std::vector<std::pair<int, int>> screen_dimensions;
        for (const TH2D* screen : measured_data) {
            screen_dimensions.emplace_back(screen->GetNbinsX(), screen->GetNbinsY());
        }
        std::vector<TH2D*> ones_projections = ForwardProject(screen_dimensions);

        // Backward project all "ones" projections to compute normalization
        for (TH2D* proj : ones_projections) {
            TH2D* ones_error = new TH2D(
                Form("ones_error_%s", proj->GetName()),
                Form("Ones Error (Screen %s)", proj->GetName()),
                proj->GetNbinsX(), proj->GetXaxis()->GetXmin(), proj->GetXaxis()->GetXmax(),
                proj->GetNbinsY(), proj->GetYaxis()->GetXmin(), proj->GetYaxis()->GetXmax()
            );
            for (int px = 0; px < ones_error->GetNbinsX(); ++px) {
                for (int py = 0; py < ones_error->GetNbinsY(); ++py) {
                    ones_error->SetBinContent(px+1, py+1, 1.0);
                }
            }
            BackwardProject({ones_error}, normalization);
            delete ones_error;
        }

        // Cleanup ones_projections
        for (TH2D* proj : ones_projections) {
            delete proj;
        }

        // SIRT iterations
        for (int iter = 0; iter < max_iterations; ++iter) {
            // Forward project current estimate to all screens
            std::vector<TH2D*> projected = ForwardProject(screen_dimensions);

            // Compute errors: measured - projected for each screen
            std::vector<TH2D*> errors;
            for (size_t screen_idx = 0; screen_idx < measured_data.size(); ++screen_idx) {
                TH2D* error = new TH2D(
                    Form("error_screen_%zu", screen_idx),
                    Form("Error (Screen %zu)", screen_idx),
                    measured_data[screen_idx]->GetNbinsX(), measured_data[screen_idx]->GetXaxis()->GetXmin(), measured_data[screen_idx]->GetXaxis()->GetXmax(),
                    measured_data[screen_idx]->GetNbinsY(), measured_data[screen_idx]->GetYaxis()->GetXmin(), measured_data[screen_idx]->GetYaxis()->GetXmax()
                );
                for (int px = 0; px < error->GetNbinsX(); ++px) {
                    for (int py = 0; py < error->GetNbinsY(); ++py) {
                        error->SetBinContent(
                            px+1, py+1,
                            measured_data[screen_idx]->GetBinContent(px+1, py+1) -
                            projected[screen_idx]->GetBinContent(px+1, py+1)
                        );
                    }
                }
                errors.push_back(error);
            }

            // Backward project all errors
            VectorVolume3D error_backprojected = reconstructed;
            for (int x = 0; x < error_backprojected.GetNx(); ++x) {
                for (int y = 0; y < error_backprojected.GetNy(); ++y) {
                    for (int z = 0; z < error_backprojected.GetNz(); ++z) {
                        error_backprojected.GetVoxel(x, y, z).value = 0.0;
                    }
                }
            }
            BackwardProject(errors, error_backprojected);

            // Update volume: x_{k+1} = x_k + λ * (A^T (b - A x_k)) / (A^T 1)
            for (int x = 0; x < reconstructed.GetNx(); ++x) {
                for (int y = 0; y < reconstructed.GetNy(); ++y) {
                    for (int z = 0; z < reconstructed.GetNz(); ++z) {
                        double norm_val = normalization.GetVoxel(x, y, z).value;
                        if (norm_val > 1e-10) { // Avoid division by zero
                            double update = relaxation * error_backprojected.GetVoxel(x, y, z).value / norm_val;
                            reconstructed.GetVoxel(x, y, z).value += update;

                            // TODO: Optionally, update the 3D vector components here.
                            // Example: reconstructed.GetVoxel(x, y, z).vector[0] += update * some_factor;
                        }
                    }
                }
            }

            // Cleanup
            for (TH2D* proj : projected) {
                delete proj;
            }
            for (TH2D* err : errors) {
                delete err;
            }

            std::cout << "Iteration " << iter+1 << "/" << max_iterations << std::endl;
        }

        return reconstructed;
    }

    // ----------------------------
    // Save the reconstructed volume to a ROOT file
    // ----------------------------
    void SaveToROOTFile(const VectorVolume3D& volume, const TString& filename) {
        TFile* file = new TFile(filename, "RECREATE");
        if (!file->IsOpen()) {
            std::cerr << "Error: Failed to open file: " << filename << std::endl;
            delete file;
            return;
        }

        // Save scalar part as TH3D
        TH3D* scalar_hist = volume.ToTH3D("scalar_volume", "Scalar Volume");
        scalar_hist->Write();

        // TODO: Optionally, save the 3D vector components as separate TH3D histograms.
        // Example:
        // TH3D* vector_x_hist = new TH3D("vector_x", "Vector X Component", ...);
        // for (int x = 0; x < volume.GetNx(); ++x) {
        //     for (int y = 0; y < volume.GetNy(); ++y) {
        //         for (int z = 0; z < volume.GetNz(); ++z) {
        //             vector_x_hist->SetBinContent(x+1, y+1, z+1, volume.GetVoxel(x, y, z).vector[0]);
        //         }
        //     }
        // }
        // vector_x_hist->Write();

        file->Close();
        delete file;
        delete scalar_hist;
        std::cout << "Volume saved to: " << filename << std::endl;
    }

private:
    VectorVolume3D fVolume; // 3D volume with 3D vectors per voxel
};
