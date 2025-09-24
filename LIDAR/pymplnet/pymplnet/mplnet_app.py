import mplnetpytools as mpt
from PyQt6.QtCore import Qt, QDir
from PyQt6.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QRadioButton,
    QListWidget,
    QWidget,
    QStyleFactory,
    QTextEdit,
    QGridLayout
)
from PyQt6.QtGui import QFont
import mplnetpytools as mpt

# Only needed for access to command line arguments
import sys

class MainWindow(QMainWindow):

    # Restart code
    EXIT_CODE_REBOOT = -12342

    def __init__(self):
        super().__init__()
        self.setWindowTitle("MPLNET GUI")
        self.vars = mpt.SelectionVariables()
        self.filevars = mpt.FileVariables()
        
        # Start with accessing the mplnet website and retrieving the ntml to parse
        url = 'https://mplnet.gsfc.nasa.gov/out/data/V3_partners/Appalachian_State/'
        html = mpt.get_mplnet_html(url)  
        self.vars.getVars(html.text)

        # Get current working directory
        qPath = QDir.currentPath()
 
        self.grid = QGridLayout()
        self.grid.setSpacing(10)

        # Information window
        infoText = '<p>Select the year, month, day, and file to download from the MPLNET website.</p>'\
                   '<p>Once the download is complete, you have the option to select a File Type '\
                   'and a Variable from that File Type to export to a .csv file.</p>'\
                   '<p>The default data directory is ../data/ from the path where this script is run.</p>'\
                   '<p>The current working directory is: ' + qPath + '</p>'

        # self.infoWindow = QLabel(infoText)
        self.infoWindow = QTextEdit(infoText)
        self.selectedText = QLabel("Selection output window")
        # Pass list to be selected from as a list of strings self.selectionWindow = SelectionWindow()
        self.selectionWindow = SelectionWindow()
        self.selectionButton = QPushButton("Begin Selection")

        # Add style properties to infoWindow
        self.infoWindow.setStyleSheet("color: black; background-color: whitesmoke; padding: 10px; border: 2px solid black")
        self.infoWindow.setAutoFillBackground(True)
        self.infoWindow.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.infoWindow.setReadOnly(True)

        # Add style properties to selectedText
        self.selectedText.setStyleSheet("color: black; background-color: whitesmoke; padding: 10px; border: 2px solid black")
        self.selectedText.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.selectedText.setWordWrap(True)

        # Add style properties to selectionWindow
        self.selectionWindow.setStyleSheet("color: black; background-color: whitesmoke; padding: 10px; border: 2px solid black")

        # Create section headers
        # Buttons work better than labels for this purpose
        self.infoLabel = QPushButton("MPLNet Data Tool Information")
        self.infoLabel.setEnabled(False)
        # Add style properties to infoLabel
        self.infoLabel.setStyleSheet("color: #ffcc00; background-color: black; padding: 5px; border: 2px solid black")

        self.selectionLabel = QPushButton("MPLNet Data Selection")
        self.selectionLabel.setEnabled(False)
        self.selectionLabel.setStyleSheet("color: #ffcc00; background-color: black; padding: 5px; border: 2px solid black")

        self.selTextLabel = QPushButton("Selection Information")
        self.selTextLabel.setEnabled(False)
        self.selTextLabel.setStyleSheet("color: #ffcc00; background-color: black; padding: 5px; border: 2px solid black")

        # Window setup with grid layout
        # Current grid is 6x4
        # Columns set to be constant size with setColumnStretch(column, size)
        # Setting row and column stretch allows for more consistent resizing
        self.grid.setColumnStretch(0, 1)
        self.grid.setColumnStretch(1, 1)
        self.grid.setColumnStretch(2, 1)
        self.grid.setColumnStretch(3, 1)
        self.grid.setColumnStretch(4, 1)
        self.grid.setRowStretch(0, 1)
        self.grid.setRowStretch(1, 1)
        self.grid.setRowStretch(2, 1)
        self.grid.setRowStretch(3, 1)
        self.grid.setRowStretch(4, 1)
        self.grid.setRowStretch(5, 1)
        self.grid.setRowStretch(6, 1)
        self.grid.setRowStretch(7, 1)

        # add grid with addWidget(widgetName, row, column, optional:rowSpan, optional:columnSpan)
        self.grid.addWidget(self.infoLabel, 1, 0, 1, 3)
        self.grid.addWidget(self.infoWindow, 2, 0, 1, 3)
        self.grid.addWidget(self.selectionLabel, 3, 0, 1, 3)
        self.grid.addWidget(self.selectionWindow, 4, 0, 2, 3)
        self.grid.addWidget(self.selectionButton, 6, 0, 1, 3)
        self.grid.addWidget(self.selTextLabel, 1, 3, 1, 2)
        self.grid.addWidget(self.selectedText, 2, 3, 4, 2)

        self.selectionButton.clicked.connect(self.selectClicked)
        self.selectionButton.clicked.connect(self.nextSelection)

        # Create restart button
        self.restartButton = QPushButton("Restart")
        self.restartButton.setEnabled(True)
        self.restartButton.clicked.connect(self.restart)
        self.grid.addWidget(self.restartButton, 0, 4)

        # Create download button
        self.downloadButton = QPushButton("Download")
        self.downloadButton.setEnabled(False)
        self.downloadButton.clicked.connect(self.downloadClicked)

        # Add downloadButton to grid
        self.grid.addWidget(self.downloadButton, 6, 3, 1, 2)

        #create the radio buttons to choose between minute, hour, and day averages
        self.minavgRadio = QRadioButton('Minute Average')
        self.minavgRadio.setEnabled(False)
        self.minavgRadio.toggled.connect(lambda:self.radioState(self.minavgRadio))
        self.grid.addWidget(self.minavgRadio, 7, 0, 1, 2)

        self.hravgRadio = QRadioButton('Hourly Average')
        self.hravgRadio.setEnabled(False)
        self.hravgRadio.toggled.connect(lambda:self.radioState(self.hravgRadio))
        self.grid.addWidget(self.hravgRadio, 7, 2, 1, 2)

        self.dayavgRadio = QRadioButton('Daily Average')
        self.dayavgRadio.setEnabled(False)
        self.dayavgRadio.toggled.connect(lambda:self.radioState(self.dayavgRadio))
        self.grid.addWidget(self.dayavgRadio, 7, 4, 1, 2)

        # Create variable selection button
        self.varSelectButton = QPushButton("Select File Type")
        self.varSelectButton.setEnabled(False)
        self.varSelectButton.clicked.connect(self.selectVars)
        self.varSelectButton.clicked.connect(self.nextVars)

        # transitButton is used to transition from downloading to file variable selection
        self.transitButton = QPushButton("Begin File Variable Selection")
        self.transitButton.setEnabled(False)
        self.transitButton.clicked.connect(self.transitFileVars)
        self.transitButton.clicked.connect(self.nextVars)

        # Create export button
        self.exportButton = QPushButton("Export to .csv")
        self.exportButton.setEnabled(False)
        self.exportButton.clicked.connect(self.exportClicked)

        # ------------------------ Build the main window ------------------------ #
        # Add layout to widget and set as central widget
        mainWidget = QWidget()

        # Add grey background to mainWidget
        mainWidget.setStyleSheet("background-color: dimgrey")
        mainWidget.setAutoFillBackground(True)

        mainWidget.setLayout(self.grid)

        # Window set to 1024x768 pixels
        mainWidget.resize(1024, 768)
        self.setCentralWidget(mainWidget)

    # ------------------------ Functions ------------------------ #
    # Function to update the selectionWindow line by line when downloading 
    # multiple files

    # Restart function
    def restart(self):
        QApplication.exit(MainWindow.EXIT_CODE_REBOOT)

    def updateSelection(self, text):
        # Get current text from selectionWindow
        currentText = []
        currentText.append(text)
        newtext = '\n'.join(currentText)
        self.selectionWindow.addItem(newtext)

    # Used to transition from downloading to file variable selection
    def transitFileVars(self):
        # Overwrite selection to allow for file variable selection
        self.grid.addWidget(self.varSelectButton, 6, 0, 1, 3)

        # Prep the selectionWindow/varSelectButton for file variable selection
        self.varSelectButton.setEnabled(True)
        self.varSelectButton.setText("Select File Type")
        self.selectionWindow.setSelectionMode(QListWidget.SelectionMode.SingleSelection)
        self.transitButton.setEnabled(False)

    def transitExport(self):
        self.grid.addWidget(self.exportButton, 6, 3, 1, 2)
        self.exportButton.setEnabled(True)
        self.minavgRadio.setEnabled(True)
        self.minavgRadio.setChecked(True)
        self.hravgRadio.setEnabled(True)
        self.dayavgRadio.setEnabled(True)

    def exportClicked(self):
        # Export the selected variable to a .csv file

        # Get the full path for files of the selected type
        if self.minavgRadio.isChecked():
            _, dirs, files = self.vars.prepDownload()
            fullpathfiles = [x + y for x, y in zip(dirs, files) if  y.find(self.filevars.selectedFileType) > -1]

            variable = self.filevars.selectedFileVars
            filename = mpt.create_export_name(self.vars, variable)

            mpt.export(filename, fullpathfiles, variable)

            self.selectionWindow.addItem("Exported " + filename)
            self.exportButton.setEnabled(False)
        if self.hravgRadio.isChecked():
            _, dirs, files = self.vars.prepDownload()
            fullpathfiles = [x + y for x, y in zip(dirs, files) if  y.find(self.filevars.selectedFileType) > -1]

            variable = self.filevars.selectedFileVars
            filename = 'HRAVG_' + mpt.create_export_name(self.vars, variable)

            mpt.HrAvg(filename, fullpathfiles, variable).to_csv(filename, header=True)

            self.selectionWindow.addItem("Exported " + filename)
            self.exportButton.setEnabled(False)
        if self.dayavgRadio.isChecked():
            _, dirs, files = self.vars.prepDownload()
            fullpathfiles = [x + y for x, y in zip(dirs, files) if  y.find(self.filevars.selectedFileType) > -1]

            variable = self.filevars.selectedFileVars
            filename = 'DAYAVG_' + mpt.create_export_name(self.vars, variable)

            mpt.DayAvg(filename, fullpathfiles, variable).to_csv(filename, header=True)

            self.selectionWindow.addItem("Exported " + filename)
            self.exportButton.setEnabled(False)

    def downloadClicked(self):
        # Download the selected files
        
        # Clear selectionWindow for downloading status
        self.selectionWindow.clear()

        # Build the download urls and update the selectionWindow
        # TODO: Add a progress bar or update files as completed for better user experience
        urls, dirs, files = self.vars.prepDownload()
        for (url, dir, file) in zip(urls, dirs, files):
            print(url +', '+ dir + file)
            text = self.vars.download(url, dir, file)
            self.updateSelection(text)

        self.downloadButton.setEnabled(False)

        # Set the file types from the download for the variable selection
        self.filevars.setFileTypes(self.vars.selectedFileTypes)

        self.grid.addWidget(self.transitButton, 6, 0, 1, 3)
        self.transitButton.setEnabled(True)

    def radioState(self,b):
        # Clear selectionWindow for downloading status
        #self.selectionWindow.clear()
        if b.text() == 'Minute Average':
            self.hravgRadio.setChecked(False)
            self.dayavgRadio.setChecked(False)
        if b.text() == 'Hourly Average':
            self.minavgRadio.setChecked(False)
            self.dayavgRadio.setChecked(False)
        if b.text() == 'Daily Average':
            self.minavgRadio.setChecked(False)
            self.hravgRadio.setChecked(False)
        

    def updateVarInfo(self):
        # Update display of variable information

        # Clear the infoWindow
        self.infoWindow.clear()

        # Show the current variable information
        self.infoWindow.setText("File Variable Information")

    def selectVars(self):
        # User select variables
        if self.filevars.peakNext():
            # get user selection and store
            varText = self.selectionWindow.selectedItems()
            if len(varText) != 0:
                # Store the selected file type
                self.filevars.storeCurrent(varText[0].text())

                # Get the file + dirs containing the variables
                _, dirs, files = self.vars.prepDownload()

                # Find the file that matches the user selection
                # Ensures that variable selection is only from the selected file type
                file = [x + y for x, y in zip(dirs, files) if  y.find(self.filevars.selectedFileType) > -1][0]

                self.filevars.setFileVars(file)

            else:
                # Handle no selection
                self.selectedText.setText("\n\n  ERROR:\n  **********************************\
                                          \nPlease select one file type\n")
                # TODO Reset selection to allow for new selection
        else: 
            self.varSelectButton.setText("Select Variable")
            # get user selection and store
            var = self.selectionWindow.selectedItems()[0].text()
            if len(var) != 0:
                # Store the selected variable
                self.filevars.storeCurrent(var)

                # Display variable description
                self.selectedText.setText(self.vars.printSelected() + self.filevars.printSelected(var))

                # Disable the varSelectButton
                self.varSelectButton.setEnabled(False)
                self.varSelectButton.setText('End variable selection')

                # Transition to export
                self.transitExport()

            else:
                # Handle no selection
                self.selectedText.setText("\n\n  ERROR:\n  **********************************\
                                          \nPlease select one variable\n")
                # TODO Reset selection to allow for new selection

    def nextVars(self):
        # Add filevars next to selectionWindow
        self.selectionWindow.clear()
        if self.filevars.peakNext():
            self.selectionWindow.addItems(self.filevars.next())

    def nextSelection(self):
        # Add vars next to selectionWindow
        self.selectionWindow.clear()
        if self.vars.peakNext():
            self.selectionWindow.addItems(self.vars.next())

    def selectClicked(self):
        # Changes button text and pulls selected items
        if self.vars.peakNext():
            self.selectionButton.setText("Next Selection")
            
            # List to hold user selected items
            listText = [item.text() for item in self.selectionWindow.selectedItems()]

            # store list into self.vars
            self.vars.storeCurrent(listText)

            # Display the entire list of selected items after storing them
            self.selectedText.setText(self.vars.printSelected())
        else:
            # When no next item exist
            # List to hold user selected items
            listText = [item.text() for item in self.selectionWindow.selectedItems()]

            # store list into self.vars
            self.vars.storeCurrent(listText)

            # Display the entire list of selected items after storing them
            self.selectedText.setText(self.vars.printSelected())

            # Ensure a selection is made in each category
            if self.vars.checkSelection():
                self.selectionButton.setText("End Selection")
                self.selectionButton.setEnabled(False)
                self.downloadButton.setEnabled(True)
            else:
                # Start selection over
                self.vars.reset()
                self.selectedText.setText("\n\n  ERROR:\n  **********************************\
                                          \nPlease select at least one item from each category\n")


class SelectionWindow(QListWidget):
    def __init__(self):
        super().__init__()

        # Set selection mode so that multiple items can be selected
        self.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)

# You need one (and only one) QApplication instance per application.
# Pass in sys.argv to allow command line arguments for your app.
# If you know you won't use command line arguments QApplication([]) works too.

# Handle restarting application
if __name__ == '__main__':
    currentExitCode = MainWindow.EXIT_CODE_REBOOT
    while currentExitCode == MainWindow.EXIT_CODE_REBOOT:
        app = QApplication([])

        # Set style to Fusion
        app.setStyle(QStyleFactory.create('Fusion'))
        app.setFont(QFont('Helvetica', 14))

        window = MainWindow()
        window.show()  # IMPORTANT!!!!! Windows are hidden by default.

        currentExitCode = app.exec()
        app = None

        # Restart commented this out
        # Start the event loop.
        # sys.exit(app.exec())

# Your application won't reach here until you exit and the event
# loop has stopped.