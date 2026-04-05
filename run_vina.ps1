# ===============================================
# Vina Docking Workflow - Single Receptor, Multi-Ligand
# Top 3 interactions only, separate cells for residue and distance
# ===============================================

# -----------------------------
# 1. User Parameters
# -----------------------------
$ligandFolder = "C:\vina\ligands"        # Ligand folder
$receptorFile = "C:\vina\receptor.pdbqt" # Single receptor in main vina folder
$outputFolder = "C:\vina\results"        # Folder to save results

$center_x = 27.555
$center_y = 25.026
$center_z = 0.111
$size_x = 40
$size_y = 40
$size_z = 40

$exhaustiveness = 8
$num_modes = 9
$distanceCutoff = 5.0     # Å cutoff for interactions
$topInteractions = 3       # Number of top closest interactions to record

$vinaExe = "C:\vina\vina.exe"
$vinaSplitExe = "C:\vina\vina_split.exe"
$csvFile = Join-Path $outputFolder "Vina_Summary_TopInteractions.csv"

# -----------------------------
# 2. Prepare output folder & CSV
# -----------------------------
if (-not (Test-Path $outputFolder)) { New-Item -ItemType Directory -Path $outputFolder }

# Initialize CSV header
$header = "Ligand,PoseFile,Energy,RMSD_Lower,RMSD_Upper"
for ($i=1; $i -le $topInteractions; $i++) {
    $header += ",Inter${i}_Res,Inter${i}_Dist"
}
$header | Out-File -FilePath $csvFile -Encoding UTF8

# -----------------------------
# 3. Function to extract top interactions
# -----------------------------
function Get-TopInteractions {
    param($ligandPoseFile, $receptorFile, $cutoff, $topN)

    $ligandAtoms = @()
    $receptorAtoms = @()
    $interactionList = @()

    # Parse ligand atoms
    foreach ($line in Get-Content $ligandPoseFile) {
        if ($line -match "^HETATM|^ATOM") {
            $ligandAtoms += [PSCustomObject]@{
                ResName = $line.Substring(17,3).Trim()
                ResNum = $line.Substring(22,4).Trim()
                X = [double]$line.Substring(30,8)
                Y = [double]$line.Substring(38,8)
                Z = [double]$line.Substring(46,8)
            }
        }
    }

    # Parse receptor atoms
    foreach ($line in Get-Content $receptorFile) {
        if ($line -match "^ATOM") {
            $receptorAtoms += [PSCustomObject]@{
                ResName = $line.Substring(17,3).Trim()
                ResNum = $line.Substring(22,4).Trim()
                X = [double]$line.Substring(30,8)
                Y = [double]$line.Substring(38,8)
                Z = [double]$line.Substring(46,8)
            }
        }
    }

    # Compute all distances
    foreach ($l in $ligandAtoms) {
        foreach ($r in $receptorAtoms) {
            $dx = $l.X - $r.X
            $dy = $l.Y - $r.Y
            $dz = $l.Z - $r.Z
            $dist = [math]::Sqrt($dx*$dx + $dy*$dy + $dz*$dz)
            if ($dist -le $cutoff) {
                $interactionList += [PSCustomObject]@{
                    ResName = $r.ResName + $r.ResNum
                    Distance = [math]::Round($dist,2)
                }
            }
        }
    }

    # Return top N closest interactions
    return ($interactionList | Sort-Object Distance | Select-Object -First $topN)
}

# -----------------------------
# 4. Docking and pose processing
# -----------------------------
$ligands = Get-ChildItem -Path $ligandFolder -Filter "*.pdbqt"

foreach ($ligand in $ligands) {
    Write-Host "`nDocking ligand:" $ligand.Name

    $ligandBase = $ligand.BaseName
    $outFile = Join-Path $outputFolder ("$ligandBase`_out.pdbqt")
    $logFile = Join-Path $outputFolder ("$ligandBase`_log.txt")
    $poseFolder = Join-Path $outputFolder ("$ligandBase`_poses")
    if (-not (Test-Path $poseFolder)) { New-Item -ItemType Directory -Path $poseFolder }

    # -----------------------------
    # Run docking
    # -----------------------------
    & $vinaExe `
        --receptor $receptorFile `
        --ligand $ligand.FullName `
        --center_x $center_x `
        --center_y $center_y `
        --center_z $center_z `
        --size_x $size_x `
        --size_y $size_y `
        --size_z $size_z `
        --out $outFile `
        --log $logFile `
        --exhaustiveness $exhaustiveness `
        --num_modes $num_modes

    if (-not (Test-Path $outFile) -or (Get-Content $outFile).Length -eq 0) {
        Write-Host "WARNING: Empty output for $ligand.Name"
        continue
    }

    # -----------------------------
    # Split poses
    # -----------------------------
    & $vinaSplitExe --input $outFile --ligand (Join-Path $poseFolder ("$ligandBase`_pose"))

    # -----------------------------
    # Parse energies, RMSD, filter poses
    # -----------------------------
    $poseFiles = Get-ChildItem -Path $poseFolder -Filter "*.pdbqt" | Sort-Object Name

    foreach ($pose in $poseFiles) {
        $lines = Get-Content $pose.FullName
        foreach ($line in $lines) {
            if ($line -match "^REMARK VINA RESULT:\s*(-?\d+\.?\d*)\s+(-?\d+\.?\d*)\s+(-?\d+\.?\d*)") {
                $energy = [double]$matches[1]
                $rmsd_lower = [double]$matches[2]
                $rmsd_upper = [double]$matches[3]

                # Skip zero or positive energy poses
                if ($energy -ge 0) { continue }

                # Extract top 3 interactions
                $topInts = Get-TopInteractions $pose.FullName $receptorFile $distanceCutoff $topInteractions

                # Build CSV row
                $row = @($ligandBase, $pose.FullName, $energy, $rmsd_lower, $rmsd_upper)
                foreach ($int in $topInts) {
                    $row += $int.ResName
                    $row += $int.Distance
                }

                # If less than 3 interactions, pad with empty cells
                $missing = $topInteractions - $topInts.Count
                for ($i=0; $i -lt $missing; $i++) {
                    $row += ""
                    $row += ""
                }

                ($row -join ",") | Out-File -FilePath $csvFile -Append -Encoding UTF8
            }
        }
    }

    Write-Host "Ligand $ligand.Name: CSV rows saved for negative-energy poses"
}

Write-Host "`nDocking complete. Summary CSV saved to $csvFile"