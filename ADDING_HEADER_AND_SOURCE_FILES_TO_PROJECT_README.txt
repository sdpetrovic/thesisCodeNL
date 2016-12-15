In this README file I will try to explain how to use header and source files (that are not classes) within a project. 


These are the steps:

1. Write a header (.h) and a source (.cpp) file.

2. Put both files in the projectLibraries folder (this is where all "libraries" are located that have been produced for this project).

3. Open the CMakeLists.txt file that is associated with this project.

4. Go to the bottom of the file and above # Add applications there should be three extra headers called 

- # Add source files.
- # Add header files.
- # Add static libraries.

5. The new files will have to be added to the source files and header files respectively. Under # Add source files
it should say set(PROJECTLIBRARIES_SOURCES
	"${SRCROOT}/projectLibraries/basicRecurrenceRelations.cpp"

	)   as an example. Simply copy the line containing the .cpp file and change the name to the new .cpp file. The same has to be done for the header file below that.

6. Re-build the CMakeLists.txt file, and now the .h and .cpp file that was added should show up in the projectLibraries folder in the master tree for the project.

7. Please make sure that the target_link_libraries part of the # Add application. includes the library project_libraries!!
