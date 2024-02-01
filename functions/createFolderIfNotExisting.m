function checkOrCreateFolder = createFolderIfNotExisting(folderpath)
    checkOrCreateFolder = isfolder(folderpath) || mkdir(folderpath);
end