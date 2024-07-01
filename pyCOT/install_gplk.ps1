# Set proper execution policy
Set-ExecutionPolicy RemoteSigned -Scope CurrentUser

# Define the URL of the file to download
$sourceUrl = ""

# Define the temporary location to store the downloaded file and the
extraction folder
$tempLocation = "$env:TEMP\glpk-4.65.zip"
$extractedFolder = "$env:TEMP\glpk-4.65"

# Define the destination system folder
$systemFolder = "${env:SystemRoot}\System32"

# Download the file
Invoke-WebRequest -Uri $sourceUrl -OutFile $tempLocation

# Extract the contents of the zip file
Expand-Archive -Path $tempLocation -DestinationPath $extractedFolder

# Copy the files from the w64 folder to the system folder
$sourceFiles = Get-ChildItem -Path "$extractedFolder\w64" -Recurse
$sourceFiles | ForEach-Object {
    $destinationPath = Join-Path $systemFolder
$_.FullName.Substring($extractedFolder.Length + 4)
    Write-Output $destinationPath
        # Copy-Item -Path $_.FullName -Destination $destinationPath
-Force
    Copy-Item -Path $_.FullName -Destination $systemFolder -Force
    
    }
    
# Clean up - delete the temporary files and folder
Remove-Item $tempLocation -Force
Remove-Item $extractedFolder -Recurse -Force
    
Write-Host "glpk-4.65 files copied to the system folder."
    
