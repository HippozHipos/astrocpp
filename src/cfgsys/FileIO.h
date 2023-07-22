#ifndef FILEIO_H
#define FILEIO_H

#include "../../util/includes.h"
namespace fs = std::filesystem;
using std::vector;
using std::string;
using std::ifstream;

class FileIO 
{
public:

    explicit FileIO(const std::string& filename)
        : fileStream(filename, std::ios::in | std::ios::out | std::ios::binary)
    {
        if (!fileStream.is_open()) 
        {
            throw std::runtime_error("Failed to open file: " + filename);
        }
    }

    ~FileIO() 
    {
        fileStream.close();
    }

    string content() 
    {
        string content;
        fileStream.seekg(0, std::ios::end);
        content.resize(fileStream.tellg());
        fileStream.seekg(0, std::ios::beg);
        fileStream.read(&content[0], content.size());
        return content;
    }

    void read()
    {
        
    }

    void write(const std::string& content) 
    {
        fileStream.seekp(0, std::ios::beg);
        fileStream.write(content.c_str(), content.size());
    }

    void File::update(const std::string& content) 
    {
        // Get the size of the file
        std::ifstream inputFile(filename, std::ios::ate | std::ios::binary);
        std::streamsize fileSize = inputFile.tellg();
        inputFile.close();

        if (fileSize > 100 * 1024 * 1024) 
        { // 100MB threshold
            // Use an alternative strategy for large files
            updateLargeFile(content);
        }
        else 
        {
            // Use the original strategy for smaller files
            updateSmallFile(content);
        }
    }

    void overwrite(const std::string& newContent) 
    {
        fileStream.close();
        fileStream.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
        if (!fileStream.is_open()) 
        {
            throw std::runtime_error("Failed to overwrite file: " + filename);
        }
        fileStream.write(newContent.c_str(), newContent.size());
    }

    void insert(const std::string& newContent, std::streampos position) 
    {
        std::string currentContent = read();
        currentContent.insert(position, newContent);
        overwrite(currentContent);
    }

private:

    void File::updateSmallFile(const std::string& content)
    {
        // 
        std::ifstream inputFile(filename);
        std::ostringstream buffer;
        buffer << inputFile.rdbuf();
        std::string fileContent = buffer.str();
        inputFile.close();

        size_t insertPos = 10; // Insert at position 10
        fileContent.insert(insertPos, content);

        std::ofstream outputFile(filename);
        outputFile << fileContent;
        outputFile.close();
    }

    void File::updateLargeFile(const std::string& content) 
    {
        // Generate a temporary file path
        std::string tempFilePath = filename + ".tmp";

        // Write the updated content to the temporary file
        std::ofstream tempFile(tempFilePath);
        tempFile << content;

        // Open the input file for reading in binary mode
        std::ifstream inputFile(filename, std::ios::binary);

        // Open the output file for writing in binary mode
        std::ofstream outputFile(filename + ".tmp", std::ios::binary);

        // Write the original file content to the output file
        outputFile << inputFile.rdbuf();

        // Close the input and output files
        inputFile.close();
        outputFile.close();

        // Rename the temporary file to replace the original file
        fs::rename(tempFilePath, filename);
    }

    std::fstream fileStream;
    const std::string filename;
    size_t linecount;
    size_t filesize;
};

#endif //FILEIO_H