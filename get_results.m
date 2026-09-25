% Boring fetch convenience script made by Gemini.
% Fetch input files from: https://zenodo.org/records/22938686

recordID = '22938686';
targetFile = 'results.zip';

% 1. Query the Zenodo REST API for the record metadata
apiUrl = sprintf('https://zenodo.org/api/records/%s', recordID);
options = weboptions('ContentType', 'json', 'Timeout', 30);
data = webread(apiUrl, options);

% 2. Extract the direct download URL for inputs.zip
downloadUrl = '';
for k = 1:numel(data.files)
    if strcmp(data.files(k).key, targetFile)
        downloadUrl = data.files(k).links.self;
        break;
    end
end

if isempty(downloadUrl)
    error('File %s not found in record %s.', targetFile, recordID);
end

% 3. Download the archive
fprintf('Downloading %s...\n', targetFile);
downloadOptions = weboptions('Timeout', 120);
savedFile = websave(targetFile, downloadUrl, downloadOptions);
fprintf('Successfully downloaded to: %s\n', savedFile);

% Optional: Unpack the archive directly
unzip(savedFile, '.');
paths;



