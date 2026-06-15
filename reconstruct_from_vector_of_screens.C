// =============================================
// SIRT Reconstruction Skeleton (ROOT C++)
// Voxel-to-Multi-Screen Forward Projection
// =============================================

#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <iostream>
#include <vector>

// ----------------------------
// Custom Forward Projection:
// Maps 3D voxels to multiple 2D screens (pixels).
// Input: 3D volume (TH3D)
// Output: Vector of 2D projections (one per screen)
// ----------------------------
std::vector<TH2D*> ForwardProject(
    const TH3D* volume,
    const std::vector<std::pair<int, int>>& screen_dimensions // (n_pixels_x, n_pixels_y) for each screen
) {
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
        for (int x = 0; x < volume->GetNbinsX(); ++x) {
            for (int y = 0; y < volume->GetNbinsY(); ++y) {
                for (int z = 0; z < volume->GetNbinsZ(); ++z) {
                    double voxel_value = volume->GetBinContent(x+1, y+1, z+1);

                    // TODO: Replace with your custom voxel-to-pixel mapping for this screen.
                    // Example: Use the 3D vector in the voxel and screen geometry to determine pixel contributions.
                    // This is a placeholder for your physics/geometry model.

                    // Example: Project voxel center to pixel coordinates for this screen.
                    // Adjust this based on your screen's position/orientation.
                    double pixel_x = (x + 0.5) * (n_pixels_x / volume->GetNbinsX());
                    double pixel_y = (y + 0.5) * (n_pixels_y / volume->GetNbinsY());

                    // Clamp to pixel bounds
                    int px = TMath::Min(TMath::Max(1, (int)(pixel_x + 0.5)), n_pixels_x);
                    int py = TMath::Min(TMath::Max(1, (int)(pixel_y + 0.5)), n_pixels_y);

                    // Accumulate voxel contribution to the pixel
                    projection->AddBinContent(px, py, voxel_value);
                }
            }
        }
    }
    return projections;
}

// ----------------------------
// Backward Projection (Adjoint):
// Maps errors from all screens back to 3D voxels.
// Input: Vector of 2D errors (one per screen)
// Output: Updates 3D volume (TH3D)
// ----------------------------
void BackwardProject(
    const std::vector<TH2D*>& errors,
    TH3D* volume
) {
    // Loop over all screens
    for (size_t screen_idx = 0; screen_idx < errors.size(); ++screen_idx) {
        const TH2D* error = errors[screen_idx];

        // Loop over all pixels in this screen
        for (int px = 0; px < error->GetNbinsX(); ++px) {
            for (int py = 0; py < error->GetNbinsY(); ++py) {
                double pixel_error = error->GetBinContent(px+1, py+1);

                // TODO: Replace with your custom pixel-to-voxel mapping for this screen.
                // Example: Distribute the pixel error to all voxels that contributed to it.
                // This is the adjoint of your forward projection for this screen.

                // Example: Map pixel back to voxels (simplified)
                for (int x = 0; x < volume->GetNbinsX(); ++x) {
                    for (int y = 0; y < volume->GetNbinsY(); ++y) {
                        for (int z = 0; z < volume->GetNbinsZ(); ++z) {
                            // TODO: Replace with your actual adjoint logic for this screen.
                            // Example: Weight the contribution based on your forward model.
                            double weight = 1.0 / (volume->GetNbinsX() * volume->GetNbinsY() * volume->GetNbinsZ());
                            volume->AddBinContent(x+1, y+1, z+1, pixel_error * weight);
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
TH3D* SIRTReconstruct(
    const TH3D* initial_volume,                     // Initial guess (3D volume)
    const std::vector<TH2D*>& measured_data,        // Measured pixel data (one TH2D per screen)
    int max_iterations = 10,
    double relaxation = 0.1                         // Relaxation parameter (0 < λ ≤ 1)
) {
    // Clone initial volume
    TH3D* reconstructed = (TH3D*)initial_volume->Clone("reconstructed");
    reconstructed->Reset();

    // Precompute normalization (sum of forward projection weights for each voxel across all screens)
    TH3D* normalization = (TH3D*)initial_volume->Clone("normalization");
    normalization->Reset();

    // Create a "ones" volume for normalization
    TH3D* ones = (TH3D*)initial_volume->Clone("ones");
    ones->Reset();
    for (int x = 0; x < ones->GetNbinsX(); ++x) {
        for (int y = 0; y < ones->GetNbinsY(); ++y) {
            for (int z = 0; z < ones->GetNbinsZ(); ++z) {
                ones->SetBinContent(x+1, y+1, z+1, 1.0);
            }
        }
    }

    // Forward project "ones" to all screens to get normalization weights
    std::vector<std::pair<int, int>> screen_dimensions;
    for (const TH2D* screen : measured_data) {
        screen_dimensions.emplace_back(screen->GetNbinsX(), screen->GetNbinsY());
    }
    std::vector<TH2D*> ones_projections = ForwardProject(ones, screen_dimensions);

    // Backward project all "ones" projections to compute normalization
    for (TH2D* proj : ones_projections) {
        TH2D* ones_error = (TH2D*)proj->Clone("ones_error");
        ones_error->Reset();
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
        std::vector<TH2D*> projected = ForwardProject(reconstructed, screen_dimensions);

        // Compute errors: measured - projected for each screen
        std::vector<TH2D*> errors;
        for (size_t screen_idx = 0; screen_idx < measured_data.size(); ++screen_idx) {
            TH2D* error = (TH2D*)measured_data[screen_idx]->Clone(Form("error_screen_%zu", screen_idx));
            error->Add(projected[screen_idx], -1.0);
            errors.push_back(error);
        }

        // Backward project all errors
        TH3D* error_backprojected = (TH3D*)initial_volume->Clone("error_bp");
        error_backprojected->Reset();
        BackwardProject(errors, error_backprojected);

        // Update volume: x_{k+1} = x_k + λ * (A^T (b - A x_k)) / (A^T 1)
        for (int x = 0; x < reconstructed->GetNbinsX(); ++x) {
            for (int y = 0; y < reconstructed->GetNbinsY(); ++y) {
                for (int z = 0; z < reconstructed->GetNbinsZ(); ++z) {
                    double norm_val = normalization->GetBinContent(x+1, y+1, z+1);
                    if (norm_val > 1e-10) { // Avoid division by zero
                        double update = relaxation * error_backprojected->GetBinContent(x+1, y+1, z+1) / norm_val;
                        reconstructed->SetBinContent(
                            x+1, y+1, z+1,
                            reconstructed->GetBinContent(x+1, y+1, z+1) + update
                        );
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
        delete error_backprojected;

        std::cout << "Iteration " << iter+1 << "/" << max_iterations << std::endl;
    }

    // Cleanup
    delete ones;
    delete normalization;

    return reconstructed;
}

// ----------------------------
// Example Usage
// ----------------------------
void RunSIRTReconstruction() {
    // 1. Create initial volume (3D)
    TH3D* initial_volume = new TH3D(
        "volume", "Initial Volume",
        32, -10, 10,  // X bins
        32, -10, 10,  // Y bins
        32, -10, 10   // Z bins
    );
    // Fill with zeros or prior
    for (int x = 0; x < 32; ++x) {
        for (int y = 0; y < 32; ++y) {
            for (int z = 0; z < 32; ++z) {
                initial_volume->SetBinContent(x+1, y+1, z+1, 0.0);
            }
        }
    }

    // 2. Create measured pixel data for multiple screens (e.g., 3 screens)
    std::vector<TH2D*> measured_data;
    for (int screen_idx = 0; screen_idx < 3; ++screen_idx) {
        TH2D* screen_data = new TH2D(
            Form("measured_screen_%d", screen_idx),
            Form("Measured Pixel Data (Screen %d)", screen_idx),
            64, 0, 64,  // X pixels
            64, 0, 64   // Y pixels
        );
        // TODO: Fill with your measured pixel data for this screen.
        measured_data.push_back(screen_data);
    }

    // 3. Run SIRT
    TH3D* result = SIRTReconstruct(initial_volume, measured_data, 10, 0.1);

    // 4. Save results
    TFile* output_file = new TFile("sirt_multi_screen_reconstruction.root", "RECREATE");
    result->Write();
    for (TH2D* screen : measured_data) {
        screen->Write();
    }
    output_file->Close();

    std::cout << "Reconstruction saved to sirt_multi_screen_reconstruction.root" << std::endl;

    // Cleanup
    delete initial_volume;
    for (TH2D* screen : measured_data) {
        delete screen;
    }
}
