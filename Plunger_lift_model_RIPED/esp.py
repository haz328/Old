import sys
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import pandas as pd
from PyQt5 import QtCore, QtGui, QtWidgets
from espSimulator import Ui_MainWindow
from expdata_input import Ui_dlg_input_expdata
from geometry_input import Ui_dlg_input_geometry
from progressBar import Ui_progressDialog
from aboutDialog import Ui_aboutDialog
from emulsionModule import Ui_dlg_emulsion_module
from erosionModule import Ui_dlg_erosion_module
import webbrowser
import espSimple
import ESPERO


# set auto scaling off
if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)


################################
# the main window of the program
################################
class MainProgram(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainProgram, self).__init__(parent)
        self.setupUi(self)
        self.move(500, 300)

        # Define variables
        self.df_catalog = pd.read_excel('Catalog.xlsx', sheetname='FLEX31')
        self.inputValues = {"R1": 0.01894, "R2": 0.039403, "TB": 0.0018875, "TV": 0.0030065,
                                "YI1": 0.015341, "YI2": 0.01046, "VOI": 0.000010506, "VOD": 0.000010108,
                                "ASF": 0, "ASB": 0, "AB": 0.0020663, "AV": 0.0025879,
                                "ADF": 0, "ADB": 0, "LI": 0.04252, "LD": 0.060315,
                                "RLK": 0.033179, "LG": 0.005, "SL": 0.000254, "EA": 0.0003,
                                "ZI": 6, "ZD": 8, "B1": 21.99, "B2": 57.03,
                                "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                                "VISW": 0.001, "ST": 0.073, "N": 3600, "SGM": 0.3,
                                "QL": 4000, "QG": 50}
        self.inputValuesCopy = {}                               # copy geometry data, dictionary
        self.df_catalogCopy = pd.DataFrame()                    # copy catalog data, pandas dataframe
        self.df_sglw = pd.DataFrame()                           # dataframe for single-phase water flow
        self.df_sglv = pd.DataFrame()                           # dataframe for single-phase viscous fluid
        self.df_flowpa = pd.DataFrame()                         # dataframe for flow pattern map
        self.df_surging = pd.DataFrame()                        # dataframe for surging
        self.df_mapping = pd.DataFrame()                        # dataframe for mapping
        self.viscosity = []                                     # viscosity list

        self.QL = np.arange(500)[1:] * 50.0                     # liquid flow rate for single-phase and mapping

        # best match flow rate input
        self.le_match.setValidator(QtGui.QDoubleValidator())
        self.le_match.setText('5500')
        self.QBEM = float(self.le_match.text())

        # open TUALP and UTULSA website
        self.label_tualp.mousePressEvent = self.label_tualp_clicked
        self.label_Utulsa.mousePressEvent = self.label_Utulsa_clicked

        # plot1 initialize
        self.figure1 = plt.figure()
        self.canvas1 = FigureCanvas(self.figure1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        self.toolbar1.setMaximumHeight(20)
        self.ax1 = self.figure1.add_subplot(111)
        plot_layout1 = QtWidgets.QVBoxLayout()
        plot_layout1.addWidget(self.toolbar1)
        plot_layout1.addWidget(self.canvas1)
        self.widget_plot1.setLayout(plot_layout1)

        # plot2 initialize
        self.figure2 = plt.figure()
        self.canvas2 = FigureCanvas(self.figure2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        self.toolbar2.setMaximumHeight(20)
        self.ax2 = self.figure2.add_subplot(111)
        plot_layout2 = QtWidgets.QVBoxLayout()
        plot_layout2.addWidget(self.toolbar2)
        plot_layout2.addWidget(self.canvas2)
        self.widget_plot2.setLayout(plot_layout2)

        # add items to combobox
        self.comboBox.addItems(["Pump1", "Pump2", "Pump3", "Pump4", "User Customize"])
        self.comboBox.setCurrentIndex(3)
        self.comboBox.currentIndexChanged.connect(self.combobox_index_changed)
        self.df_catalogCopy = self.df_catalog
        self.df_catalogCopy.columns = ['N', 'Q', 'P']

        self.btn_read_file.clicked.connect(self.btn_read_file_clicked)
        self.btn_geom_edit.clicked.connect(self.btn_geom_edit_clicked)
        self.btn_geom_apply.clicked.connect(self.btn_geom_apply_clicked)

        # water performance group box
        self.btn_water_edit.clicked.connect(self.btn_water_edit_clicked)
        self.btn_water_view.clicked.connect(self.btn_water_view_clicked)

        # Tuning group box
        self.btn_clear_tune.clicked.connect(self.btn_clear_tune_clicked)
        self.btn_tune.clicked.connect(self.btn_tune_clicked)

        # Single-phase group box
        self.btn_clear_single_phase.clicked.connect(self.btn_clear_single_phase_clicked)
        self.btn_add_viscosity.clicked.connect(self.btn_add_viscosity_clicked)
        self.btn_single_phase_calculate.clicked.connect(self.btn_single_phase_calculate_clicked)
        self.btn_emulsion.clicked.connect(self.btn_emulsion_clicked)
        self.btn_erosion.clicked.connect(self.btn_erosion_clicked)

        # Two-phase group box
        self.btn_surging_calc.clicked.connect(self.btn_surging_calc_clicked)
        self.btn_mapping_calc.clicked.connect(self.btn_mapping_calc_clicked)
        self.btn_flowpa.clicked.connect(self.btn_flowpa_clicked)
        self.btn_clear_two_phase.clicked.connect(self.btn_clear_two_phase_clicked)

        #  Data Output group box
        self.radioButton_single_phase.toggled.connect(self.radioButton_single_phase_toggled)
        self.radioButton_surging.toggled.connect(self.radioButton_surging_toggled)
        self.radioButton_mapping.toggled.connect(self.radioButton_mapping_toggled)
        self.radioButton_flow_pattern.toggled.connect(self.radioButton_flow_pattern_toggled)
        self.btn_out_clear.clicked.connect(self.btn_out_clear_clicked)
        self.btn_out_2_figure.clicked.connect(self.btn_out_2_figure_clicked)
        self.btn_out_2_file.clicked.connect(self.btn_out_2_file_clicked)
        self.tableWidget.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

        # save and exit buttons
        self.btn_exit.clicked.connect(self.btn_exit_clicked)
        self.btn_save.clicked.connect(self.btn_save_clicked)

        # Menu-File
        self.actionOpen.triggered.connect(self.btn_read_file_clicked)
        self.actionSave.triggered.connect(self.btn_save_clicked)
        self.actionExit.triggered.connect(QtWidgets.qApp.quit)

        # Menu-Data
        self.actionImport_Geomerty.triggered.connect(self.btn_read_file_clicked)
        self.actionImport_Water_Curve.triggered.connect(self.btn_water_edit_clicked)
        self.actionExport_Data.triggered.connect(self.btn_out_2_file_clicked)
        self.actionExport_Figure.triggered.connect(self.btn_out_2_figure_clicked)
        self.actionAdd_Viscosity.triggered.connect(self.btn_add_viscosity_clicked)

        # Menu-Run
        self.actionTune.triggered.connect(self.btn_tune_clicked)
        self.actionSingle_phase.triggered.connect(self.btn_single_phase_calculate_clicked)
        self.actionSurging.triggered.connect(self.btn_surging_calc_clicked)
        self.actionMapping.triggered.connect(self.btn_mapping_calc_clicked)
        self.actionFlow_Pattern_Prediction.triggered.connect(self.btn_flowpa_clicked)
        self.actionEmulsion_Module.triggered.connect(self.btn_emulsion_clicked)
        self.actionErosion_Module.triggered.connect(self.btn_erosion_clicked)

        # Menu-Help
        self.actionAbout.triggered.connect(self.actionAbout_triggered)
        self.actionHelp.triggered.connect(self.actionHelp_triggered)

    #########################
    # open web
    @QtCore.pyqtSlot()
    def label_tualp_clicked(self, e):
        webbrowser.open("http://www.tualp.utulsa.edu/")

    @QtCore.pyqtSlot()
    def label_Utulsa_clicked(self, e):
        webbrowser.open("https://utulsa.edu/")

    #########################
    # INPUTS group box slots
    # pump model selection -> combobox on the main window
    @QtCore.pyqtSlot()
    def combobox_index_changed(self):
        if self.comboBox.currentIndex() == 1:
            # DN1750-Pump2
            self.inputValues = {"R1": 1.9875E-2, "R2": 3.5599E-2, "TB": 1.7E-3, "TV": 3.12E-3,
                                "YI1": 1.3536E-2, "YI2": 7.13E-3, "VOI": 6.283E-6, "VOD": 7.063E-6,
                                "ASF": 6.8159E-04, "ASB": 6.549E-04, "AB": 6.9356E-04, "AV": 7.1277E-04,
                                "ADF": 1.0605E-03, "ADB": 6.6436E-04, "LI": 3.9E-02, "LD": 5.185E-02,
                                "RLK": 0.04, "LG": 0.01, "SL": 0.00005, "EA": 0.0002,
                                "ZI": 6, "ZD": 8, "B1": 20.3, "B2": 36.2,
                                "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                                "VISW": 0.001, "ST": 0.073, "N": 3500, "SGM": 0.3,
                                "QL": 1750, "QG": 50}
            self.df_catalog = pd.read_excel('Catalog.xlsx', sheetname='DN1750')

        # TE2700-Pump1
        elif self.comboBox.currentIndex() == 0:
            self.inputValues = {"R1": 0.017496, "R2": 0.056054, "TB": 0.00272, "TV": 0.00448,
                                "YI1": 0.012194, "YI2": 0.007835, "VOI": 0.000016119, "VOD": 0.000011153,
                                "ASF": 0.004659, "ASB": 0.003718, "AB": 0.001319, "AV": 0.001516,
                                "ADF": 0.000832, "ADB": 0.00137, "LI": 0.076, "LD": 0.08708,
                                "RLK": 0.056209, "LG": 0.00806, "SL": 0.00005, "EA": 0.0002,
                                "ZI": 5, "ZD": 9, "B1": 29.5, "B2": 34.7,
                                "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                                "VISW": 0.001, "ST": 0.073, "N": 3500, "SGM": 0.3,
                                "QL": 2700, "QG": 50}
            self.df_catalog = pd.read_excel('Catalog.xlsx', sheetname='TE2700')

        # GC6100-Pump3
        elif self.comboBox.currentIndex() == 2:
            self.inputValues = {"R1": 0.014351, "R2": 0.050013, "TB": 0.0025896, "TV": 0.002894,
                                "YI1": 0.017399, "YI2": 0.013716, "VOI": 1.512E-5, "VOD": 1.9818E-5,
                                "ASF": 8.9654E-4, "ASB": 9.4143E-4, "AB": 1.0333E-3, "AV": 1.769E-3,
                                "ADF": 2.0486E-3, "ADB": 1.0301E-3, "LI": 0.0529, "LD": 0.0839,
                                "RLK": 6.1237E-2, "LG": 0.0015475, "SL": 3.81E-4, "EA": 0.0003,
                                "ZI": 7, "ZD": 8, "B1": 33.375, "B2": 41.387,
                                "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                                "VISW": 0.001, "ST": 0.073, "N": 3600, "SGM": 0.3,
                                "QL": 6100, "QG": 50}
            self.df_catalog = pd.read_excel('Catalog.xlsx', sheetname='GC6100')

        # FLEX31-pump4
        elif self.comboBox.currentIndex() == 3:
            self.inputValues = {"R1": 0.01894, "R2": 0.039403, "TB": 0.0018875, "TV": 0.0030065,
                                "YI1": 0.015341, "YI2": 0.01046, "VOI": 0.000010506, "VOD": 0.000010108,
                                "ASF": 0, "ASB": 0, "AB": 0.0020663, "AV": 0.0025879,
                                "ADF": 0, "ADB": 0, "LI": 0.04252, "LD": 0.060315,
                                "RLK": 0.033179, "LG": 0.005, "SL": 0.000254, "EA": 0.0003,
                                "ZI": 6, "ZD": 8, "B1": 21.99, "B2": 57.03,
                                "DENL": 1000, "DENG": 11.2, "VISL": 0.001, "VISG": 0.000018,
                                "VISW": 0.001, "ST": 0.073, "N": 3600, "SGM": 0.3,
                                "QL": 4000, "QG": 50}
            self.df_catalog = pd.read_excel('Catalog.xlsx', sheetname='FLEX31')

        # user customized pump geometry
        else:
            self.inputValues = {"R1": 0,    "R2": 0,    "TB": 0,    "TV": 0,    "YI1": 0,   "YI2": 0,   "VOI": 0,
                                "VOD": 0,   "ASF": 0,   "ASB": 0,   "AB": 0,    "AV": 0,    "ADF": 0,   "ADB": 0,
                                "LI": 0,    "LD": 0,    "RLK": 0,   "LG": 0,    "SL": 0,    "EA": 0,    "ZI": 0,
                                "ZD": 0,    "B1": 0,    "B2": 0,    "DENL": 0,  "DENG": 0,  "VISL": 0,  "VISG": 0,
                                "VISW": 0,  "ST": 0,    "N": 0,     "SGM": 0,   "QL": 0,   "QG": 0}
            self.df_catalog = pd.DataFrame()
            dlg = DlgInputGeomData(self.inputValues)
            if dlg.exec_():
                self.inputValues = dlg.inputValues.copy()

        try:
            self.df_catalogCopy = self.df_catalog
            self.df_catalogCopy.columns = ['N', 'Q', 'P']
        except ValueError:
            pass

    # read in geometry file -> read file button on the main window
    @QtCore.pyqtSlot()
    def btn_read_file_clicked(self):
        if self.comboBox.currentIndex() != 3:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Pump Model is selected!")
            msg.setWindowTitle("Error")
            msg.exec_()
        else:
            options = QtWidgets.QFileDialog.Options()
            options |= QtWidgets.QFileDialog.DontUseNativeDialog
            filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Read Input File", "",
                                                                "All Files (*);;txt Files (*.txt);;csv Files (*.csv)",
                                                                options=options)
            if filename:
                if filename.split('.')[-1] == 'txt':
                    df_input = pd.read_table(filename).drop(["mark"], 1).set_index(["item"])
                    for index in df_input.index:
                        self.inputValues[index] = df_input.loc[index].value.item()
                elif filename.split('.')[-1] == 'csv':
                    df_input = pd.read_csv(filename).drop(["mark"], 1).set_index(["item"])
                    for index in df_input.index:
                        self.inputValues[index] = df_input.loc[index].value.item()
                else:
                    msg = QtWidgets.QMessageBox()
                    msg.setIcon(QtWidgets.QMessageBox.Critical)
                    msg.setText("Please read the correct data file!")
                    msg.setWindowTitle("File Open Error")
                    msg.exec_()

    # edit input parameters
    @QtCore.pyqtSlot()
    def btn_geom_edit_clicked(self):
        dlg = DlgInputGeomData(self.inputValues)
        if dlg.exec_():
            self.inputValues = dlg.inputValues.copy()

    # apply the changes made in edit
    @QtCore.pyqtSlot()
    def btn_geom_apply_clicked(self):
        self.inputValuesCopy = self.inputValues.copy()
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("Input parameters are done!")
        msg.setWindowTitle("Inputs check")
        msg.exec_()

    # input experimental data for single-phase water performance
    @QtCore.pyqtSlot()
    def btn_water_edit_clicked(self):
        dlg = DlgInputExpData()
        dlg.table_input.setColumnCount(3)
        dlg.table_input.setRowCount(len(self.df_catalog) + 1)

        # fill the default values
        for i in range(len(self.df_catalog.index)):
            for j in range(len(self.df_catalog.columns)):
                dlg.table_input.setItem(i, j, QtWidgets.QTableWidgetItem(str(self.df_catalog.iloc[i, j])))

        # retrieve all data from QTable
        if dlg.exec_():
            self.df_catalog = pd.DataFrame([], index=range(dlg.table_input.rowCount()),
                                           columns=range(dlg.table_input.columnCount()))

            for i in range(dlg.table_input.rowCount()):
                for j in range(dlg.table_input.columnCount()):
                    try:
                        self.df_catalog.iloc[i, j] = float(dlg.table_input.item(i, j).text())
                    except:
                        pass

            self.df_catalog.columns = ['N', 'Q', 'P']

            if self.df_catalog['N'].unique().size > 2:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Please use unique rotational speed!")
                msg.setWindowTitle("Error")
                msg.exec_()

            self.df_catalogCopy = self.df_catalog.drop(['N'], 1).dropna()

            if self.df_catalogCopy.shape[0] < 5:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Warning)
                msg.setText("Please add more experimental data!")
                msg.setWindowTitle("Warning")
                msg.exec_()

    # plot single-phase water performance
    @QtCore.pyqtSlot()
    def btn_water_view_clicked(self):
        self.figure1.clear()
        self.ax1 = self.figure1.add_subplot(111)
        self.ax1.tick_params(direction='in', labelsize=6)
        self.ax1.set_title('Single-phase Performance', fontsize=6)
        self.ax1.set_xlabel(r'$Q_L$ (bpd)', fontsize=6)
        self.ax1.set_ylabel('P (psi)', fontsize=6)
        self.ax1.plot(self.df_catalogCopy['Q'], self.df_catalogCopy['P'], 'bo', markerfacecolor="None", label='Catalog')
        self.ax1.legend(frameon=False, fontsize=6)
        self.canvas1.draw()

    #############################
    # Calculation group box slots
    # Tune
    @QtCore.pyqtSlot()
    def btn_clear_tune_clicked(self):
        self.figure1.clear()
        self.canvas1.draw()

    @QtCore.pyqtSlot()
    def btn_tune_clicked(self):
        index = 0
        hpsgl = []
        hesgl = []
        heesgl = []
        hfsgl = []
        htsgl = []
        hresgl = []
        qlksgl = []
        self.df_sglw = pd.DataFrame()
        self.QBEM = float(self.le_match.text())
        inputs = self.inputValues.copy()
        inputs['VISL'] = inputs['VISW']         # Always use water for tuning!

        self.ax1 = self.figure1.add_subplot(111)
        self.ax1.plot(self.df_catalogCopy['Q'], self.df_catalogCopy['P'], 'bo', markerfacecolor="None", label='Catalog')
        self.ax1.tick_params(direction='in', labelsize=6)
        self.ax1.set_title('Single-phase Performance', fontsize=6)
        self.ax1.set_xlabel('Q (bpd)', fontsize=6)
        self.ax1.set_ylabel('P (psi)', fontsize=6)

        for q in self.QL:
            HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = self.single_phase_calculation(inputs,
                                                                                                      self.QBEM, q)
            if HP_sgl > 0:
                index += 1
                hpsgl.append(HP_sgl)
                hesgl.append(HE_sgl)
                heesgl.append(HEE_sgl)
                hfsgl.append(HF_sgl)
                htsgl.append(HT_sgl)
                hresgl.append(HRE_sgl)
                qlksgl.append(QLK_sgl)
            else:
                break
        # save the single-phase water flow dataframe
        self.df_sglw = pd.DataFrame({'QL': self.QL[:len(hpsgl)], 'HP': hpsgl, 'HE': hesgl, 'HEE': heesgl, 'HF': hfsgl,
                                     'HT': htsgl, 'HRE': hresgl, 'QLK': qlksgl})

        if self.cb_hold_plots.isChecked():
            self.ax1.plot(self.df_sglw.QL, self.df_sglw.HP, 'b-')
            self.canvas1.draw()

        else:
            self.figure1.clear()
            ax = self.figure1.add_subplot(111)
            ax.plot(self.df_catalogCopy['Q'], self.df_catalogCopy['P'], 'bo', markerfacecolor="None", label='Catalog')
            ax.plot(self.df_sglw.QL,self.df_sglw.HP, 'b-', label='$Q_{BEM}=$'+'{}bpd'.format(int(self.QBEM)))
            ax.tick_params(direction='in', labelsize=6)
            ax.set_title('Single-phase Performance', fontsize=6)
            ax.set_xlabel(r'$Q_L$ (bpd)', fontsize=6)
            ax.set_ylabel('P (psi)', fontsize=6)
            ax.legend(frameon=False, fontsize=6)
            self.canvas1.draw()

        # fill the data output table
        headers = ['QL (bpd)', 'HP (psi)', 'HE (psi)', 'HEE (psi)', 'HF (psi)', 'HT (psi)',
                   'HRE (psi)', 'QLK (bpd)']
        columns = ['QL', 'HP', 'HE', 'HEE', 'HF', 'HT', 'HRE', 'QLK']
        self.fill_table(self.tableWidget, self.df_sglw, headers, columns)

    # Single-phase
    @QtCore.pyqtSlot()
    def btn_clear_single_phase_clicked(self):
        self.viscosity = []
        self.figure1.clear()
        self.canvas1.draw()

    @QtCore.pyqtSlot()
    def btn_add_viscosity_clicked(self):
        # add input variables
        if np.abs(self.inputValues['VISW'] - self.inputValues['VISL']) > 1:
            if self.inputValues['VISW'] not in self.viscosity:
                self.viscosity.append(self.inputValues['VISW'])

        if len(self.viscosity) < 5:
            d, okPressed = QtWidgets.QInputDialog.getDouble(self, "Input viscosity", "Value (cp):", 1, 0, 10000, 1)
            if okPressed and d not in self.viscosity:
                self.viscosity.append(d)
        else:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setText("More than 5 viscosities have been added!")
            msg.setWindowTitle("Input check")
            msg.exec_()

    @QtCore.pyqtSlot()
    def btn_single_phase_calculate_clicked(self):
        viscosity = np.sort(self.viscosity)
        self.QBEM = float(self.le_match.text())
        if len(viscosity) < 1:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Information)
            msg.setText("Please add more viscosity!")
            msg.setWindowTitle("Input check")
            msg.exec_()
        else:
            if self.df_sglw.shape[0] < 5:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("Best match flow rate has not been tuned!")
                msg.setWindowTitle("Calculation check")
                msg.exec_()
            else:
                self.figure1.clear()
                ax = self.figure1.add_subplot(111)
                ax.plot(self.df_catalogCopy['Q'], self.df_catalogCopy['P'], 'bo',
                        markerfacecolor="None", label='Catalog')
                ax.plot(self.df_sglw.QL, self.df_sglw.HP, 'b-',
                        label='$Q_{BEM}=$' + '{}bpd'.format(int(self.QBEM)))
                index = 0
                for vi in viscosity:
                    hpsgl = []
                    hesgl = []
                    heesgl = []
                    hfsgl = []
                    htsgl = []
                    hresgl = []
                    qlksgl = []
                    self.df_sglv = pd.DataFrame()
                    inputs = self.inputValues.copy()
                    inputs['VISL'] = vi/1000.0        #convert cp to pa.s
                    colors = ['r', 'g', 'k', 'y', 'm']
                    for q in self.QL:
                        HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = self.single_phase_calculation(
                            inputs, self.QBEM, q)
                        if HP_sgl > 0:
                            hpsgl.append(HP_sgl)
                            hesgl.append(HE_sgl)
                            heesgl.append(HEE_sgl)
                            hfsgl.append(HF_sgl)
                            htsgl.append(HT_sgl)
                            hresgl.append(HRE_sgl)
                            qlksgl.append(QLK_sgl)
                        else:
                            break
                    # save the single-phase water flow dataframe
                    self.df_sglv = pd.DataFrame(
                        {'QL': self.QL[:len(hpsgl)], 'HP': hpsgl, 'HE': hesgl, 'HEE': heesgl, 'HF': hfsgl,
                        'HT': htsgl, 'HRE': hresgl, 'QLK': qlksgl})

                    ax.plot(self.df_sglv.QL, self.df_sglv.HP, c=colors[index], label='$\mu$={}cp'.format(vi))
                    index += 1

                ax.tick_params(direction='in', labelsize=6)
                ax.set_title('Single-phase Performance', fontsize=6)
                ax.set_xlabel(r'$Q_L$ (bpd)', fontsize=6)
                ax.set_ylabel('P (psi)', fontsize=6)
                ax.legend(frameon=False, fontsize=6)
                self.canvas1.draw()

    @QtCore.pyqtSlot()
    def btn_emulsion_clicked(self):
        self.QBEM = float(self.le_match.text())
        dlg = DlgEmulsionModule(self.inputValues.copy(), self.df_catalogCopy, self.QBEM)
        if dlg.exec_():
            pass

    @QtCore.pyqtSlot()
    def btn_erosion_clicked(self):
        self.QBEM = float(self.le_match.text())
        dlg = DlgErosionModule(self.inputValues.copy(), self.df_catalogCopy, self.QBEM)
        if dlg.exec_():
            pass

    # Two-phase
    @QtCore.pyqtSlot()
    def btn_surging_calc_clicked(self):
        self.df_surging = pd.DataFrame()
        ql = self.inputValues['QL']
        self.QBEM = float(self.le_match.text())

        HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = self.single_phase_calculation(self.inputValues,
                                                                                                  self.QBEM, ql)
        if self.df_sglw.shape[0] < 5:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Best match flow rate has not been tuned!")
            msg.setWindowTitle("Calculation check")
            msg.exec_()
        else:
            if HP_sgl < 0:
                msg = QtWidgets.QMessageBox()
                msg.setIcon(QtWidgets.QMessageBox.Critical)
                msg.setText("The liquid flow rate is unphysical or the pump performance is negative")
                msg.setWindowTitle("Calculation check")
                msg.exec_()
            else:
                progbar = SurgingProgressBar(self.inputValues, self.QBEM)
                if progbar.exec_():
                    self.df_surging = progbar.df

            self.figure2.clear()
            ax = self.figure2.add_subplot(111)
            ax.plot(self.df_surging['GVF'], self.df_surging['P_single-phase'], 'k-', label='Single-phase')
            ax.plot(self.df_surging['GVF'], self.df_surging['P_homogeneous'], 'b--', label='Homogeneous')
            ax.plot(self.df_surging['GVF'], self.df_surging['P_model'], 'r-.', label='Mechanistic')
            ax.tick_params(direction='in', labelsize=6)
            ax.set_title(r'Surging performance $Q_L$ = {} bpd'.format(ql), fontsize=6)
            ax.set_xlabel(r'$\lambda_G$ (%)', fontsize=6)
            ax.set_ylabel('P (psi)', fontsize=6)

            axplus = ax.twinx()
            axplus.plot(self.df_surging['GVF'], self.df_surging['GVF'], 'g:', label='Homogeneous')
            axplus.plot(self.df_surging['GVF'], self.df_surging['GV'],  'm:', label='Mechanistic')

            axplus.set_ylabel(r'$\alpha_G$ (%)', fontsize=6)
            axplus.tick_params(direction='in', labelsize=6)

            ax.legend(frameon=False, loc='center left', title='Left Y', fontsize=6)
            axplus.legend(frameon=False, loc='center right', title='Right Y', fontsize=6)
            self.canvas2.draw()

            # output table
            headers = ['GVF (%)', 'GV (%)', 'P_single-phase (psi)', 'P_homogeneous (psi)', 'P_model (psi)', 'Flow_pattern']
            columns = ['GVF', 'GV', 'P_single-phase', 'P_homogeneous', 'P_model', 'Flow_pattern']
            self.fill_table(self.tableWidget, self.df_surging, headers, columns)

    @QtCore.pyqtSlot()
    def btn_mapping_calc_clicked(self):
        self.df_mapping = pd.DataFrame()
        self.QBEM = float(self.le_match.text())
        qg = self.inputValues['QG']

        if self.df_sglw.shape[0] < 5:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Best match flow rate has not been tuned!")
            msg.setWindowTitle("Calculation check")
            msg.exec_()
        else:
            progbar = MappingProgressBar(self.inputValues, self.QBEM)
            if progbar.exec_():
                self.df_mapping = progbar.df

            self.figure2.clear()
            ax = self.figure2.add_subplot(111)
            ax.plot(self.df_mapping['QL'], self.df_mapping['P_single-phase'], 'b-', label='Single-phase')
            ax.plot(self.df_mapping['QL'], self.df_mapping['P_homogeneous'], 'g--', label='Homogeneous')
            ax.plot(self.df_mapping['QL'], self.df_mapping['P_model'], 'r-.', label='Mechanistic')

            axplus = ax.twinx()
            axplus.plot(self.df_mapping['QL'], self.df_mapping['GVF'], 'g:', label='Homogeneous')
            axplus.plot(self.df_mapping['QL'], self.df_mapping['GV'], 'm:', label='Mechanistic')

            ax.tick_params(direction='in', labelsize=6)
            ax.set_title(r'Mapping Performance $Q_G$ = {} bpd'.format(qg), fontsize=6)
            ax.set_xlabel(r'$Q_L$ (bpd)', fontsize=6)
            ax.set_ylabel('P (psi)', fontsize=6)
            ax.legend(frameon=False, title='Left Y', loc='center', fontsize=6)

            axplus.tick_params(direction='in', labelsize=6)
            axplus.legend(frameon=False, loc='upper right', title='Right Y', fontsize=6)
            axplus.set_ylabel(r'$\alpha_G$ (%)', fontsize=6)

            self.canvas2.draw()

            # output table
            headers = ['QL (bpd)', 'GVF (%)', 'GV (%)', 'P_single-phase (psi)', 'P_homogeneous (psi)', 'P_model (psi)',
                       'Flow_pattern']
            columns = ['QL', 'GVF', 'GV', 'P_single-phase', 'P_homogeneous', 'P_model', 'Flow_pattern']
            self.fill_table(self.tableWidget, self.df_mapping, headers, columns)

    @QtCore.pyqtSlot()
    def btn_flowpa_clicked(self):
        self.df_flowpa = pd.DataFrame()
        self.QBEM = float(self.le_match.text())
        self.QBEM = float(self.le_match.text())

        if self.df_sglw.shape[0] < 5:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Best match flow rate has not been tuned!")
            msg.setWindowTitle("Calculation check")
            msg.exec_()
        else:
            progbar = FlowPaProgressBar(self.inputValues, self.QBEM)
            if progbar.exec_():
                self.df_flowpa = progbar.df

            self.figure2.clear()
            ax = self.figure2.add_subplot(111)
            b2, = ax.plot(self.df_flowpa['QSG2'], self.df_flowpa['QL'], 's', markerfacecolor='none', label='BUB to INT')
            b3, = ax.plot(self.df_flowpa['QSG3'], self.df_flowpa['QL'], 'd', markerfacecolor='none', label='INT to SEG')
            b1, = ax.plot(self.df_flowpa['QSG1'], self.df_flowpa['QL'], 'o', markerfacecolor='none', label='DB to BUB')
            ax.set_xlim(0, np.square(self.inputValues['N']/3200)*1000)
            ax.set_xlabel(r'$Q_{G}$ (bpd)', fontsize=6)
            ax.set_ylabel(r'$Q_{L}$ (bpd)', fontsize=6)
            ax.tick_params(direction='in', labelsize=6)
            ax.legend([b1, b2, b3], ['DB to BUB', 'BUB to INT', 'INT to SEG'], frameon=False, fontsize=6)
            ax.set_title('Flow pattern Map', fontsize=6)
            self.canvas2.draw()

            # fill output table
            headers = ['QL (bpd)', 'QSG1 (bpd)', 'QSG2 (bpd)', 'QSG3 (bpd)']
            columns = ['QL', 'QSG1', 'QSG2', 'QSG3']
            self.fill_table(self.tableWidget, self.df_flowpa, headers, columns)

    @QtCore.pyqtSlot()
    def btn_clear_two_phase_clicked(self):
        self.figure2.clear()
        self.canvas2.draw()


    #############################
    # Data output group box slots
    @QtCore.pyqtSlot()
    def radioButton_single_phase_toggled(self):
        if self.radioButton_single_phase.isChecked():
            if self.df_sglw.shape[0] < 5:
                pass
            else:
                headers = ['QL (bpd)', 'HP (psi)', 'HE (psi)', 'HEE (psi)', 'HF (psi)', 'HT (psi)',
                           'HRE (psi)', 'QLK (bpd)']
                columns = ['QL', 'HP', 'HE', 'HEE', 'HF', 'HT', 'HRE', 'QLK']
                self.fill_table(self.tableWidget, self.df_sglw, headers, columns)
                self.figure1.clear()
                ax = self.figure1.add_subplot(111)
                ax.plot(self.df_catalogCopy['Q'], self.df_catalogCopy['P'], 'bo', markerfacecolor="None", label='Catalog')
                ax.plot(self.df_sglw.QL, self.df_sglw.HP, 'b-', label='$Q_{BEM}=$' + '{}bpd'.format(int(self.QBEM)))
                ax.tick_params(direction='in', labelsize=6)
                ax.set_title('Single-phase Performance', fontsize=6)
                ax.set_xlabel(r'$Q_L$ (bpd)', fontsize=6)
                ax.set_ylabel('P (psi)', fontsize=6)
                ax.legend(frameon=False, fontsize=6)
                self.canvas1.draw()

    @QtCore.pyqtSlot()
    def radioButton_surging_toggled(self):
        if self.radioButton_surging.isChecked():
            if self.df_surging.shape[0] < 5:
                pass
            else:
                headers = ['GVF (%)', 'GV (%)', 'P_single-phase (psi)', 'P_homogeneous (psi)', 'P_model (psi)',
                           'Flow_pattern']
                columns = ['GVF', 'GV', 'P_single-phase', 'P_homogeneous', 'P_model', 'Flow_pattern']
                self.fill_table(self.tableWidget, self.df_surging, headers, columns)
                self.figure2.clear()
                ax = self.figure2.add_subplot(111)
                ax.plot(self.df_surging['GVF'], self.df_surging['P_single-phase'], 'k-', label='Single-phase')
                ax.plot(self.df_surging['GVF'], self.df_surging['P_homogeneous'], 'b--', label='Homogeneous')
                ax.plot(self.df_surging['GVF'], self.df_surging['P_model'], 'r-.', label='Mechanistic')
                ax.tick_params(direction='in', labelsize=6)
                ax.set_title(r'Surging performance $Q_L$ = {} bpd'.format(self.inputValues['QL']), fontsize=6)
                ax.set_xlabel(r'$\lambda_G$ (%)', fontsize=6)
                ax.set_ylabel('P (psi)', fontsize=6)

                axplus = ax.twinx()
                axplus.plot(self.df_surging['GVF'], self.df_surging['GVF'], 'g:', label='Homogeneous')
                axplus.plot(self.df_surging['GVF'], self.df_surging['GV'], 'm:', label='Mechanistic')
                axplus.set_ylabel(r'$\alpha_G$ (%)', fontsize=6)
                axplus.tick_params(direction='in', labelsize=6)

                ax.legend(frameon=False, loc='center left', title='Left Y', fontsize=6)
                axplus.legend(frameon=False, loc='center right', title='Right Y', fontsize=6)
                self.canvas2.draw()

    @QtCore.pyqtSlot()
    def radioButton_mapping_toggled(self):
        if self.radioButton_mapping.isChecked():
            if self.df_mapping.shape[0] < 5:
                pass
            else:
                headers = ['QL (bpd)', 'GVF (%)', 'GV (%)', 'P_single-phase (psi)', 'P_homogeneous (psi)',
                           'P_model (psi)',
                           'Flow_pattern']
                columns = ['QL', 'GVF', 'GV', 'P_single-phase', 'P_homogeneous', 'P_model', 'Flow_pattern']
                self.fill_table(self.tableWidget, self.df_mapping, headers, columns)

                self.figure2.clear()
                ax = self.figure2.add_subplot(111)
                ax.plot(self.df_mapping['QL'], self.df_mapping['P_single-phase'], 'b-', label='Single-phase')
                ax.plot(self.df_mapping['QL'], self.df_mapping['P_homogeneous'], 'g--', label='Homogeneous')
                ax.plot(self.df_mapping['QL'], self.df_mapping['P_model'], 'r-.', label='Mechanistic')

                axplus = ax.twinx()
                axplus.plot(self.df_mapping['QL'], self.df_mapping['GVF'], 'g:', label='Homogeneous')
                axplus.plot(self.df_mapping['QL'], self.df_mapping['GV'], 'm:', label='Mechanistic')

                ax.tick_params(direction='in', labelsize=6)
                ax.set_title(r'Mapping Performance $Q_G$ = {} bpd'.format(self.inputValues['QG']), fontsize=6)
                ax.set_xlabel('Q_L (bpd)', fontsize=6)
                ax.set_ylabel('P (psi)', fontsize=6)
                ax.legend(frameon=False, title='Left Y', loc='center', fontsize=6)
                axplus.tick_params(direction='in', labelsize=6)
                axplus.legend(frameon=False, loc='upper right', title='Right Y', fontsize=6)
                axplus.set_ylabel(r'$\alpha_G$ (%)', fontsize=6)
                self.canvas2.draw()

    @QtCore.pyqtSlot()
    def radioButton_flow_pattern_toggled(self):
        if self.radioButton_flow_pattern.isChecked():
            if self.df_flowpa.shape[0] < 5:
                pass
            else:
                headers = ['QL (bpd)', 'QSG1 (bpd)', 'QSG2 (bpd)', 'QSG3 (bpd)']
                columns = ['QL', 'QSG1', 'QSG2', 'QSG3']
                self.fill_table(self.tableWidget, self.df_flowpa, headers, columns)

                self.figure2.clear()
                ax = self.figure2.add_subplot(111)
                b2, = ax.plot(self.df_flowpa['QSG2'], self.df_flowpa['QL'], 's', markerfacecolor='none',
                              label='Boundary II')
                b3, = ax.plot(self.df_flowpa['QSG3'], self.df_flowpa['QL'], 'd', markerfacecolor='none',
                              label='Boundary III')
                b1, = ax.plot(self.df_flowpa['QSG1'], self.df_flowpa['QL'], 'o', markerfacecolor='none',
                              label='Boundary I')
                ax.set_xlim(0, np.square(self.inputValues['N'] / 3200) * 1000)
                ax.set_xlabel(r'$Q_{G}$ (bpd)', fontsize=6)
                ax.set_ylabel(r'$Q_{L}$ (bpd)', fontsize=6)
                ax.tick_params(direction='in', labelsize=6)
                ax.legend([b1, b2, b3], ['DB to BUB', 'BUB to INT', 'INT to SEG'], frameon=False, fontsize=6)
                ax.set_title('Flow pattern Map', fontsize=6)
                self.canvas2.draw()

    @QtCore.pyqtSlot()
    def btn_out_clear_clicked(self):
        self.tableWidget.clear()

    @QtCore.pyqtSlot()
    def btn_out_2_figure_clicked(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName1, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Figure1", "",
                                                            "All Files (*);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
                                                            "PDF Files (*.pdf);;PNG Files (*.png)", options=options)
        if fileName1:
            self.figure1.savefig(fileName1, dpi=300)

        fileName2, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Figure2", "",
                                                             "All Files (*);;BMP Files (*.bmp);;JPG Files (*.jpg);;"
                                                             "PDF Files (*.pdf);;PNG Files (*.png)", options=options)
        if fileName2:
            self.figure2.savefig(fileName2, dpi=300)

    @QtCore.pyqtSlot()
    def btn_out_2_file_clicked(self):
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save Current Data Table", "",
                                                            "All Files (*);;CSV Files (*.csv);;"
                                                            "Excel Files (*.xls);;Text Files (*.txt)", options=options)
        columns = []
        try:
            for i in range(self.tableWidget.columnCount()):
                columns.append(str(self.tableWidget.horizontalHeaderItem(i).text()))

            df = pd.DataFrame({}, index=range(self.tableWidget.rowCount()), columns=range(self.tableWidget.columnCount()))
            for r in range(self.tableWidget.rowCount()):
                for c in range(self.tableWidget.columnCount()):
                    df.iloc[r, c] = self.tableWidget.item(r, c).text()

            df.columns = columns
            if fileName:
                if str(fileName).split('.')[1] == 'xls':
                    df.to_excel(fileName)
                else:
                    df.to_csv(fileName)

        except AttributeError:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("No Data Available in Output Table!")
            msg.setWindowTitle("Error")
            msg.exec_()

    #############
    # SAVE & EXIT
    @QtCore.pyqtSlot()
    def btn_exit_clicked(self):
        QtWidgets.qApp.quit()

    @QtCore.pyqtSlot()
    def btn_save_clicked(self):
        self.figure1.savefig('Figure1.jpg', dpi=300)
        self.figure2.savefig('Figure2.jpg', dpi=300)

        writer = pd.ExcelWriter('Output.xlsx', engine='xlsxwriter')
        self.df_sglw.to_excel(writer, sheet_name='Single-phase water performance')
        self.df_surging.to_excel(writer, sheet_name='Surging performance')
        self.df_mapping.to_excel(writer, sheet_name='Mapping performance')
        self.df_flowpa.to_excel(writer, sheet_name='Flow pattern map')
        writer.save()

        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The calculation data and figures have been saved!")
        msg.setWindowTitle("Save Files")
        msg.exec_()

    # Help Menu
    @QtCore.pyqtSlot()
    def actionAbout_triggered(self):
        dlg = DlgAbout()
        if dlg.exec_():
            dlg.show()

    @QtCore.pyqtSlot()
    def actionHelp_triggered(self):
        pass

    ##########################################
    # Static methods in the Main program class
    @staticmethod
    def single_phase_calculation(pump, qbem, ql):
        """
        :param pump:pump dictionary, refer to MainProgram __init__ function
        :param qbem:best match flow rate (bpd)
        :param ql:  liquid flow rate (bpd)
        :return:    HP, HE, HEE, HF, HT, HRE, QLK (pressure in psi, flow rate in bpd)
        """
        R1 = pump['R1']
        R2 = pump['R2']
        TB = pump['TB']
        TV = pump['TV']
        YI1 = pump['YI1']
        YI2 = pump['YI2']
        VOI = pump['VOI']
        VOD = pump['VOD']
        ASF = pump['ASF']
        ASB = pump['ASB']
        AB = pump['AB']
        AV = pump['AV']
        ADF = pump['ADF']
        ADB = pump['ADB']
        LI = pump['LI']
        LD = pump['LD']
        RLK = pump['RLK']
        LG = pump['LG']
        SL = pump['SL']
        EA = pump['EA']
        ZI = pump['ZI']
        ZD = pump['ZD']
        B1 = pump['B1']
        B2 = pump['B2']
        DENL = pump['DENL']
        VISL = pump['VISL']
        VISW = pump['VISW']
        N = pump['N']
        SGM = pump['SGM']
        QBEM = qbem

        ql = ql * 0.15898 / 24.0 / 3600.0
        QBEM = QBEM * 0.15898 / 24.0 / 3600.0
        #QBEM = QBEM * 0.15898 / 24.0 / 3600.0 * (N / 3500.0) / np.power(VISL / VISW, 0.05)

        HE, HEE, HP, HIO, HFI, HFD, HTI, HTD, HLK, U1, U2, W1, W2, C1, C2, C2E, QLK, DI = \
            espSimple.sgl_mod.sgl(ql, QBEM, SGM, R1, R2, RLK, N, B1, B2, TB, TV, YI1, YI2, VOI, VOD, AB,
                                  AV, ASF, ASB, ADF, ADB, EA, LI, LD, ZI, ZD, LG, SL, DENL, VISL, VISW)

        # return pressure in psi, flow rate in bpd
        HP = HP * 9.8 * DENL / 6891.2
        HE = HE * 9.8 * DENL / 6891.2
        HEE = HEE * 9.8 * 997.0 / 6891.2
        HF = (HFI + HFD) * 9.8 * DENL / 6891.2
        HT = (HTI + HTD) * 9.8 * DENL / 6891.2
        HRE = abs(HE - HEE)
        QLK = QLK * 24.0 * 3600.0 / 0.15898
        return HP, HE, HEE, HF, HT, HRE, QLK

    @staticmethod
    def two_phase_calculation(pump, qbem, ql, qg):
        """
        :param pump:pump dictionary, refer to MainProgram __init__ function
        :param qbem:best match flow rate (bpd)
        :param ql:  liquid flow rate (bpd)
        :param qg:  gas flow rate (bpd)
        :return:    GF, GV, PP, FGL (pressure in psi)
        """
        R1 = pump['R1']
        R2 = pump['R2']
        TB = pump['TB']
        TV = pump['TV']
        YI1 = pump['YI1']
        YI2 = pump['YI2']
        VOI = pump['VOI']
        VOD = pump['VOD']
        ASF = pump['ASF']
        ASB = pump['ASB']
        AB = pump['AB']
        AV = pump['AV']
        ADF = pump['ADF']
        ADB = pump['ADB']
        LI = pump['LI']
        LD = pump['LD']
        RLK = pump['RLK']
        LG = pump['LG']
        SL = pump['SL']
        EA = pump['EA']
        ZI = pump['ZI']
        ZD = pump['ZD']
        B1 = pump['B1']
        B2 = pump['B2']
        DENL = pump['DENL']
        DENG = pump['DENG']
        VISL = pump['VISL']
        VISG = pump['VISG']
        VISW = pump['VISW']
        ST = pump['ST']
        N = pump['N']
        SGM = pump['SGM']
        QBEM = qbem

        ql = ql * 0.15898 / 24.0 / 3600.0
        QBEM = QBEM * 0.15898 / 24.0 / 3600.0
        #QBEM = QBEM * 0.15898 / 24.0 / 3600.0 * (N / 3500.0) / np.power(VISL / VISW, 0.05)
        qg = qg * 0.15898 / 24.0 / 3600.0
        GF = qg / (qg + ql)

        HE, HEE, HP, HIO, HFI, HFD, HTI, HTD, HLK, U1, U2, W1, W2, C1, C2, C2E, QLK, DI = \
            espSimple.sgl_mod.sgl(ql, QBEM, SGM, R1, R2, RLK, N, B1, B2, TB, TV, YI1, YI2, VOI, VOD, AB,
                               AV, ASF, ASB, ADF, ADB, EA, LI, LD, ZI, ZD, LG, SL, DENL, VISL, VISW)

        GV, PE, PEE, PP, PIO, PFI, PFD, PTI, PTD, PLK, U1, U2, W1, W2, C1, C2, C2E, QLK, VSR, DB, CD, FGL \
            = espSimple.gl_mod.gl(HP, ql, qg, QBEM, SGM, GF, R1, R2, RLK, N, B1, B2, TB, TV, YI1, YI2, VOI, VOD, AB, AV, ASF,
                               ASB, ADF, ADB, EA, LI, LD, LG, SL, ZI, ZD, DENL, DENG, VISL, VISW, VISG, ST)

        return GF, GV, PP, FGL

    @staticmethod
    def flow_pattern_calculation(pump, qbem, ql):
        """
        :param pump:pump dictionary, refer to MainProgram __init__ function
        :param qbem:best match flow rate (bpd)
        :param ql:  liquid flow rate (bpd)
        :return:    gvf1, gvf2, gvf3 at the transition boundary corresponding to a fix liquid flow rate (bpd)
        """
        G = 9.81
        PI = 3.14159265359
        R1 = pump['R1']
        R2 = pump['R2']
        TB = pump['TB']
        TV = pump['TV']
        YI1 = pump['YI1']
        YI2 = pump['YI2']
        VOI = pump['VOI']
        VOD = pump['VOD']
        ASF = pump['ASF']
        ASB = pump['ASB']
        AB = pump['AB']
        AV = pump['AV']
        ADF = pump['ADF']
        ADB = pump['ADB']
        LI = pump['LI']
        LD = pump['LD']
        RLK = pump['RLK']
        LG = pump['LG']
        SL = pump['SL']
        EA = pump['EA']
        ZI = pump['ZI']
        ZD = pump['ZD']
        B1 = pump['B1']
        B2 = pump['B2']
        DENL = pump['DENL']
        DENG = pump['DENG']
        VISL = pump['VISL']
        VISG = pump['VISG']
        VISW = pump['VISW']
        ST = pump['ST']
        N = pump['N']
        SGM = pump['SGM']
        QBEM = qbem

        # dPump = 2.0 * R2
        RI = (R1 + R2) / 2.0
        YI = (YI1 + YI2) / 2.0
        Omega = 2.0 * PI * N / 60.0
        AI = VOI / LI
        # AD = VOD / LD
        DI = 4.0 * VOI / (AB + ASF + ASB)
        # DD = 4.0 * VOD / (AV + ADF + ADB)
        EDI = EA / DI

        ql = ql * 0.15898 / 24.0 / 3600.0
        QBEM = QBEM * 0.15898 / 24.0 / 3600.0 * (N / 3500.0) / np.power(VISL / VISW, 0.05)

        HE, HEE, HP, HIO, HFI, HFD, HTI, HTD, HLK, U1, U2, W1, W2, C1, C2, C2E, QLK, DI = \
            espSimple.sgl_mod.sgl(ql, QBEM, SGM, R1, R2, RLK, N, B1, B2, TB, TV, YI1, YI2, VOI, VOD, AB,
                                  AV, ASF, ASB, ADF, ADB, EA, LI, LD, ZI, ZD, LG, SL, DENL, VISL, VISW)

        gvf1 = espSimple.gl_mod.get_lamdacrit1(DENG, DENL, G, HP, Omega, ql, QLK, RI, ST, VOI, ZI)
        gvf2 = espSimple.gl_mod.get_lamdacrit2(DENL, DENG, G, HP, N, Omega, PI, ql, QLK, RI, ST, VISL, VOI, YI, ZI, TB)
        gvf3 = espSimple.gl_mod.get_lamdacrit3(AI, B1, B2, DI, EDI, ql, DENL, G, Omega, LI, DENG, PI, RI, VISL, VISG,
                                               ST, gvf1, gvf2, ZI, QLK)

        return gvf1, gvf2, gvf3

    @staticmethod
    def fill_table(qtable, df_data, headers, columns):
        """
        :param qtable:  Qtable to fill
        :param df_data: dataframe
        :param headers: headers to change table headers
        :param columns: columns for selecting data
        :return:        None
        """
        qtable.clear()
        df = df_data[columns]
        qtable.setColumnCount(df.shape[1])
        qtable.setRowCount(df.shape[0])
        qtable.setHorizontalHeaderLabels(headers)
        qtable.resizeColumnsToContents()
        for i in range(df.shape[0]):
            for j in range(df.shape[1]):
                try:
                    item = str(np.around(df.iloc[i, j], decimals=3))
                except TypeError:
                    item = str(df.iloc[i, j])
                qtable.setItem(i, j, QtWidgets.QTableWidgetItem(item))


#####################################################################
# class for input geometry data dialog
#####################################################################
class DlgInputGeomData(QtWidgets.QDialog, Ui_dlg_input_geometry):
    def __init__(self, initvalues, parent=None):
        super(DlgInputGeomData, self).__init__(parent)
        self.setupUi(self)
        # initialize input blanks
        self.le_R1.setValidator(QtGui.QDoubleValidator())
        self.le_R1.setText(str(initvalues["R1"]))
        self.le_R2.setValidator(QtGui.QDoubleValidator())
        self.le_R2.setText(str(initvalues["R2"]))
        self.le_TB.setValidator(QtGui.QDoubleValidator())
        self.le_TB.setText(str(initvalues["TB"]))
        self.le_TV.setValidator(QtGui.QDoubleValidator())
        self.le_TV.setText(str(initvalues["TV"]))
        self.le_YI1.setValidator(QtGui.QDoubleValidator())
        self.le_YI1.setText(str(initvalues["YI1"]))
        self.le_YI2.setValidator(QtGui.QDoubleValidator())
        self.le_YI2.setText(str(initvalues["YI2"]))
        self.le_VOI.setValidator(QtGui.QDoubleValidator())
        self.le_VOI.setText(str(initvalues["VOI"]))
        self.le_VOD.setValidator(QtGui.QDoubleValidator())
        self.le_VOD.setText(str(initvalues["VOD"]))
        self.le_ASF.setValidator(QtGui.QDoubleValidator())
        self.le_ASF.setText(str(initvalues["ASF"]))
        self.le_ASB.setValidator(QtGui.QDoubleValidator())
        self.le_ASB.setText(str(initvalues["ASB"]))
        self.le_AB.setValidator(QtGui.QDoubleValidator())
        self.le_AB.setText(str(initvalues["AB"]))
        self.le_AV.setValidator(QtGui.QDoubleValidator())
        self.le_AV.setText(str(initvalues["AV"]))
        self.le_ADF.setValidator(QtGui.QDoubleValidator())
        self.le_ADF.setText(str(initvalues["ADF"]))
        self.le_ADB.setValidator(QtGui.QDoubleValidator())
        self.le_ADB.setText(str(initvalues["ADB"]))
        self.le_LI.setValidator(QtGui.QDoubleValidator())
        self.le_LI.setText(str(initvalues["LI"]))
        self.le_LD.setValidator(QtGui.QDoubleValidator())
        self.le_LD.setText(str(initvalues["LD"]))
        self.le_RLK.setValidator(QtGui.QDoubleValidator())
        self.le_RLK.setText(str(initvalues["RLK"]))
        self.le_LG.setValidator(QtGui.QDoubleValidator())
        self.le_LG.setText(str(initvalues["LG"]))
        self.le_SL.setValidator(QtGui.QDoubleValidator())
        self.le_SL.setText(str(initvalues["SL"]))
        self.le_EA.setValidator(QtGui.QDoubleValidator())
        self.le_EA.setText(str(initvalues["EA"]))
        self.le_ZI.setValidator(QtGui.QDoubleValidator())
        self.le_ZI.setText(str(initvalues["ZI"]))
        self.le_ZD.setValidator(QtGui.QDoubleValidator())
        self.le_ZD.setText(str(initvalues["ZD"]))
        self.le_B1.setValidator(QtGui.QDoubleValidator())
        self.le_B1.setText(str(initvalues["B1"]))
        self.le_B2.setValidator(QtGui.QDoubleValidator())
        self.le_B2.setText(str(initvalues["B2"]))
        self.le_DENL.setValidator(QtGui.QDoubleValidator())
        self.le_DENL.setText(str(initvalues["DENL"]))
        self.le_DENG.setValidator(QtGui.QDoubleValidator())
        self.le_DENG.setText(str(initvalues["DENG"]))
        self.le_VISL.setValidator(QtGui.QDoubleValidator())
        self.le_VISL.setText(str(initvalues["VISL"]))
        self.le_VISG.setValidator(QtGui.QDoubleValidator())
        self.le_VISG.setText(str(initvalues["VISG"]))
        self.le_VISW.setValidator(QtGui.QDoubleValidator())
        self.le_VISW.setText(str(initvalues["VISW"]))
        self.le_ST.setValidator(QtGui.QDoubleValidator())
        self.le_ST.setText(str(initvalues["ST"]))
        self.le_N.setValidator(QtGui.QDoubleValidator())
        self.le_N.setText(str(initvalues["N"]))
        self.le_SGM.setValidator(QtGui.QDoubleValidator())
        self.le_SGM.setText(str(initvalues["SGM"]))
        self.le_QL.setValidator(QtGui.QDoubleValidator())
        self.le_QL.setText(str(initvalues["QL"]))
        self.le_QG.setValidator(QtGui.QDoubleValidator())
        self.le_QG.setText(str(initvalues["QG"]))

        self.inputValues = initvalues.copy()
        self.buttonBox.accepted.connect(self.ok_pressed)

    @QtCore.pyqtSlot()
    def ok_pressed(self):
        try:
            self.inputValues["R1"] = float(self.le_R1.text())
            self.inputValues["R2"] = float(self.le_R2.text())
            self.inputValues["TB"] = float(self.le_TB.text())
            self.inputValues["TV"] = float(self.le_TV.text())
            self.inputValues["YI1"] = float(self.le_YI1.text())
            self.inputValues["YI2"] = float(self.le_YI2.text())
            self.inputValues["VOI"] = float(self.le_VOI.text())
            self.inputValues["VOD"] = float(self.le_VOD.text())
            self.inputValues["ASF"] = float(self.le_ASF.text())
            self.inputValues["ASB"] = float(self.le_ASB.text())
            self.inputValues["AB"] = float(self.le_AB.text())
            self.inputValues["AV"] = float(self.le_AV.text())
            self.inputValues["ADF"] = float(self.le_ADF.text())
            self.inputValues["ADB"] = float(self.le_ADB.text())
            self.inputValues["LI"] = float(self.le_LI.text())
            self.inputValues["LD"] = float(self.le_LD.text())
            self.inputValues["RLK"] = float(self.le_RLK.text())
            self.inputValues["LG"] = float(self.le_LG.text())
            self.inputValues["SL"] = float(self.le_SL.text())
            self.inputValues["EA"] = float(self.le_EA.text())
            self.inputValues["ZI"] = float(self.le_ZI.text())
            self.inputValues["ZD"] = float(self.le_ZD.text())
            self.inputValues["B1"] = float(self.le_B1.text())
            self.inputValues["B2"] = float(self.le_B2.text())
            self.inputValues["DENL"] = float(self.le_DENL.text())
            self.inputValues["DENG"] = float(self.le_DENG.text())
            self.inputValues["VISL"] = float(self.le_VISL.text())
            self.inputValues["VISG"] = float(self.le_VISG.text())
            self.inputValues["VISW"] = float(self.le_VISW.text())
            self.inputValues["ST"] = float(self.le_ST.text())
            self.inputValues["N"] = float(self.le_N.text())
            self.inputValues["SGM"] = float(self.le_SGM.text())
            self.inputValues["QL"] = float(self.le_QL.text())
            self.inputValues["QG"] = float(self.le_QG.text())
        except ValueError:
            QtWidgets.QMessageBox.warning(self, "Wrong Input", "Please Check Input values!")


####################################################################
# class for input experimental data manually dialog
####################################################################
class DlgInputExpData(QtWidgets.QDialog, Ui_dlg_input_expdata):
    def __init__(self, parent=None):
        super(DlgInputExpData, self).__init__(parent)
        self.setupUi(self)

        self.btn_accept.clicked.connect(self.accept)
        self.btn_cancel.clicked.connect(self.reject)
        self.btn_read.clicked.connect(self.btn_read_clicked)
        self.btn_add.clicked.connect(self.btn_add_clicked)
        self.btn_delete.clicked.connect(self.btn_delete_clicked)
        self.btn_clear.clicked.connect(self.btn_clear_clicked)

    @QtCore.pyqtSlot()
    def btn_read_clicked(self):
        buttonReply = QtWidgets.QMessageBox.question(self, 'Read Data File', "Do you want to overwrite current data?",
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                     QtWidgets.QMessageBox.No)
        if buttonReply == QtWidgets.QMessageBox.Yes:
            self.table_input.clearContents()
            options = QtWidgets.QFileDialog.Options()
            options |= QtWidgets.QFileDialog.DontUseNativeDialog
            filename, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Read Data File", "",
                                                                "All Files (*);;txt Files (*.txt);;csv Files (*.csv)",
                                                                options=options)
            df_data = pd.DataFrame()
            if filename:
                if filename.split('.')[-1] == 'txt':
                    df_data = pd.read_table(filename)
                elif filename.split('.')[-1] == 'csv':
                    df_data = pd.read_csv(filename)
                else:
                    msg = QtWidgets.QMessageBox()
                    msg.setIcon(QtWidgets.QMessageBox.Critical)
                    msg.setText("Please read the correct data file!")
                    msg.setWindowTitle("File Open Error")
                    msg.exec_()

            self.table_input.setRowCount(df_data.shape[0])
            self.table_input.setColumnCount(df_data.shape[1])

            for i in range(df_data.shape[0]):
                for j in range(df_data.shape[1]):
                    self.table_input.setItem(i, j, QtWidgets.QTableWidgetItem(str(df_data.iloc[i, j])))

    @QtCore.pyqtSlot()
    def btn_add_clicked(self):
        rowcount = self.table_input.rowCount()
        self.table_input.setRowCount(rowcount + 1)

    @QtCore.pyqtSlot()
    def btn_delete_clicked(self):
        rowcount = self.table_input.rowCount()
        self.table_input.setRowCount(rowcount - 1)

    @QtCore.pyqtSlot()
    def btn_clear_clicked(self):
        buttonReply = QtWidgets.QMessageBox.question(self, 'Clear Data', "Do you want to clear input data?",
                                                     QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                                                     QtWidgets.QMessageBox.No)
        if buttonReply == QtWidgets.QMessageBox.Yes:
            self.table_input.clearContents()


#########################
# class for about dialog
#########################
class DlgAbout(QtWidgets.QDialog, Ui_aboutDialog):
    def __init__(self):
        super(DlgAbout, self).__init__()
        self.setupUi(self)
        self.pushButton.clicked.connect(self.accept)


###################################
# classes for progress bar
###################################
class SurgingProgressBar(QtWidgets.QDialog, Ui_progressDialog):
    def __init__(self, inputValues, QBEM, parent=None):
        super(SurgingProgressBar, self).__init__(parent)
        self.setupUi(self)
        self.inputValues = inputValues
        self.QBEM = QBEM
        self.ql = self.inputValues['QL']
        self.df = pd.DataFrame()                    # a DataFrame to store results
        self.n = 30                                 # maximum data count in progressbar

        self.qg_max = 0.5 * self.ql                 # construct an array for surging iteration with gas flow rates
        self.QG = np.arange(self.n) * (self.qg_max / self.n)
        self.QG[0] = 5.0

        self.progressBar.setValue(0)
        self.progressBar.setTextVisible(True)
        self.surgingThread = SurgingThread(self.n, self.inputValues, self.QBEM, self.QG)
        self.surgingThread.total.connect(self.progressBar.setMaximum)
        self.surgingThread.update.connect(self.update)
        self.surgingThread.finished.connect(self.finished)
        self.btn_OK.clicked.connect(self.accept)

        self.icon = 0
        self.surgingThread.start()

    @QtCore.pyqtSlot()
    def finished(self):
        self.df = pd.DataFrame({'GVF': self.surgingThread.GF, 'GV': self.surgingThread.GV,
                                'P_single-phase': self.surgingThread.PP_sgl, 'P_homogeneous': self.surgingThread.PP_homo,
                                'P_model': self.surgingThread.PP, 'Flow_pattern': self.surgingThread.FGL})
        self.progressBar.setValue(self.n)
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The calculation is finished!")
        msg.setWindowTitle("Finish")
        msg.exec_()

    @QtCore.pyqtSlot()
    def update(self):
        self.icon += 1
        self.progressBar.setValue(self.icon)


class MappingProgressBar(QtWidgets.QDialog, Ui_progressDialog):
    def __init__(self, inputValues, QBEM, parent=None):
        super(MappingProgressBar, self).__init__(parent)
        self.setupUi(self)
        self.inputValues = inputValues
        self.QBEM = QBEM
        self.qg = self.inputValues['QG']
        self.df = pd.DataFrame()                    # a DataFrame to store results

        self.progressBar.setValue(0)
        self.progressBar.setTextVisible(True)
        self.MappingThread = MappingThread(self.inputValues, self.QBEM)
        self.MappingThread.total.connect(self.progressBar.setMaximum)
        self.MappingThread.update.connect(self.update)
        self.MappingThread.finished.connect(self.finished)
        self.btn_OK.clicked.connect(self.accept)

        self.icon = 0
        self.MappingThread.start()

    @QtCore.pyqtSlot()
    def finished(self):
        self.df = pd.DataFrame({'QL': self.MappingThread.QL, 'GVF': self.MappingThread.GF, 'GV': self.MappingThread.GV,
                                'P_single-phase': self.MappingThread.PP_sgl, 'P_homogeneous': self.MappingThread.PP_homo,
                                'P_model': self.MappingThread.PP, 'Flow_pattern': self.MappingThread.FGL})
        self.progressBar.setValue(self.MappingThread.n)
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The calculation is finished!")
        msg.setWindowTitle("Finish")
        msg.exec_()

    @QtCore.pyqtSlot()
    def update(self):
        self.icon += 1
        self.progressBar.setValue(self.icon)


class FlowPaProgressBar(QtWidgets.QDialog, Ui_progressDialog):
    def __init__(self, inputValues, QBEM, parent=None):
        super(FlowPaProgressBar, self).__init__(parent)
        self.setupUi(self)
        self.inputValues = inputValues
        self.QBEM = QBEM
        self.df = pd.DataFrame()                    # a DataFrame to store results

        self.progressBar.setValue(0)
        self.progressBar.setTextVisible(True)
        self.FlowPaThread = FlowPaThread(self.inputValues, self.QBEM)
        self.FlowPaThread.total.connect(self.progressBar.setMaximum)
        self.FlowPaThread.update.connect(self.update)
        self.FlowPaThread.finished.connect(self.finished)
        self.btn_OK.clicked.connect(self.accept)

        self.icon = 0
        self.FlowPaThread.start()

    @QtCore.pyqtSlot()
    def finished(self):
        self.df = pd.DataFrame({'QL': self.FlowPaThread.QL,     'QSG1': self.FlowPaThread.QSG1,
                                'QSG2': self.FlowPaThread.QSG2, 'QSG3': self.FlowPaThread.QSG3})
        self.progressBar.setValue(self.FlowPaThread.n)
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The calculation is finished!")
        msg.setWindowTitle("Finish")
        msg.exec_()

    @QtCore.pyqtSlot()
    def update(self):
        self.icon += 1
        self.progressBar.setValue(self.icon)


###################################
#  class for multi-thread
###################################
class SurgingThread(QtCore.QThread):
    # class variables for signals
    total = QtCore.pyqtSignal(object)
    update = QtCore.pyqtSignal()

    def __init__(self, n, inputValues, QBEM, QG):
        super(SurgingThread, self).__init__()
        self.n = n                          # maximum run times
        self.inputValues = inputValues      # input dictionary
        self.QBEM = QBEM
        self.QG = QG                        # QG array for surging calculation
        self.ql = self.inputValues['QL']
        self.PP_homo = []                   # homogeneous model for comparison
        self.PP_sgl = []                    # single-phase performance
        self.GF = []
        self.GV = []
        self.PP = []
        self.FGL = []

        HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = \
            MainProgram.single_phase_calculation(self.inputValues, self.QBEM, self.ql)

        self.HP_sgl = HP_sgl

    def run(self):
        self.total.emit(self.n)

        i = 0
        while i <= self.n:
            self.msleep(1)          # A necessary pause to show the progress
            qg = self.QG[i]
            gf, gv, pp, fgl = MainProgram.two_phase_calculation(self.inputValues, self.QBEM, self.ql, qg)
            if gf <= 0.3:   # GVF is high enough to make the plot
                self.PP_sgl.append(self.HP_sgl)
                self.GF.append(gf * 100)
                self.GV.append(gv * 100)
                self.FGL.append(fgl)
                c_h = (1 - gf) + gf * self.inputValues['DENG'] / self.inputValues['DENL']
                c_m = (1 - gv) + gv * self.inputValues['DENG'] / self.inputValues['DENL']
                self.PP_homo.append(c_h * self.HP_sgl)
                self.PP.append(c_m * self.HP_sgl)
            else:
                break
            i += 1
            self.update.emit()


class MappingThread(QtCore.QThread):
    # class variables for signals
    total = QtCore.pyqtSignal(object)
    update = QtCore.pyqtSignal()

    def __init__(self, inputValues, QBEM):
        super(MappingThread, self).__init__()
        self.inputValues = inputValues      # input dictionary
        self.QBEM = QBEM
        self.qg = self.inputValues['QG']
        self.QL, self.PP_sgl = self.get_ql()
        self.n = len(self.QL)
        self.PP_homo = []                   # homogeneous model for comparison
        self.GF = []
        self.GV = []
        self.PP = []
        self.FGL = []

    # construct a list/array of QL for iteration
    def get_ql(self):
        QL = []
        PP_sgl = []
        ql = 50
        HP_sgl = 1
        while HP_sgl > 0:
            HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = \
                MainProgram.single_phase_calculation(self.inputValues, self.QBEM, ql)
            QL.append(ql)
            PP_sgl.append(HP_sgl)
            ql += 50
        return QL[:-1], PP_sgl[:-1]

    def run(self):
        self.total.emit(self.n)
        i = 0
        while i < self.n:
            self.msleep(1)          # A necessary pause to show the progress
            ql = self.QL[i]
            gf, gv, pp, fgl = MainProgram.two_phase_calculation(self.inputValues, self.QBEM, ql, self.qg)
            if pp < 0:
                break
            self.GF.append(gf * 100)
            self.GV.append(gv * 100)
            self.PP.append(pp / 6891.2)
            self.FGL.append(fgl)
            c_h = (1 - gf) + gf * self.inputValues['DENG'] / self.inputValues['DENL']
            self.PP_homo.append(c_h * self.PP_sgl[i])
            i += 1
            self.update.emit()

        # reset the QL, PP_sgl index range to match GF, GV, PP, FGL range
        self.QL = self.QL[:i]
        self.PP_sgl = self.PP_sgl[:i]


class FlowPaThread(QtCore.QThread):
    # class variables for signals
    total = QtCore.pyqtSignal(object)
    update = QtCore.pyqtSignal()

    def __init__(self, inputValues, QBEM):
        super(FlowPaThread, self).__init__()
        self.inputValues = inputValues  # input dictionary
        self.QBEM = QBEM
        self.QL = self.get_ql()
        self.n = len(self.QL)
        self.QSG1 = []
        self.QSG2 = []
        self.QSG3 = []

    # construct a list/array of QL for iteration
    def get_ql(self):
        QL = []
        ql = 50
        HP_sgl = 1
        while HP_sgl > 0:
            HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = \
                MainProgram.single_phase_calculation(self.inputValues, self.QBEM, ql)
            QL.append(ql)
            ql += 50
        return QL[:-1]

    def run(self):
        self.total.emit(self.n)
        i = 0
        while i < self.n:
            self.msleep(1)  # A necessary pause to show the progress
            ql = self.QL[i]
            gvf1, gvf2, gvf3 = MainProgram.flow_pattern_calculation(self.inputValues, self.QBEM, ql)
            if gvf1 > gvf2:
                gvf2 = gvf1
            self.QSG1.append(gvf1 / (1.0 - gvf1) * ql)
            self.QSG2.append(gvf2 / (1.0 - gvf2) * ql)
            self.QSG3.append(gvf3 / (1.0 - gvf3) * ql)
            i += 1
            self.update.emit()


#######################################
# class for emulsion module calculation
#######################################
class DlgEmulsionModule(QtWidgets.QDialog, Ui_dlg_emulsion_module):
    def __init__(self, pump, catalog_data, QBEM, parent=None):
        super(DlgEmulsionModule, self).__init__(parent)
        self.setupUi(self)

        # initialize new dialogue inputs
        self.le_oil_viscosity.setValidator(QtGui.QDoubleValidator())
        self.le_oil_viscosity.setText(str(0.1))
        self.le_oil_density.setValidator(QtGui.QDoubleValidator())
        self.le_oil_density.setText(str(850.))
        self.le_water_viscosity.setValidator(QtGui.QDoubleValidator())
        self.le_water_viscosity.setText(str(0.001))
        self.le_water_density.setValidator(QtGui.QDoubleValidator())
        self.le_water_density.setText(str(1000.))
        self.le_water_cut.setValidator(QtGui.QDoubleValidator())
        self.le_water_cut.setText(str(15))
        self.le_surface_tension.setValidator(QtGui.QDoubleValidator())
        self.le_surface_tension.setText(str(0.02))
        self.le_N.setValidator(QtGui.QDoubleValidator())
        self.le_N.setText(str(3500.))
        self.le_Q.setValidator(QtGui.QDoubleValidator())
        self.le_Q.setText(str(1000.))
        self.le_stage_number.setValidator(QtGui.QIntValidator())
        self.le_stage_number.setText(str(1))

        self.oil_viscosity = float(self.le_oil_viscosity.text())
        self.oil_density = float(self.le_oil_density.text())
        self.water_viscosity = float(self.le_water_viscosity.text())
        self.water_density = float(self.le_water_density.text())
        self.water_cut = float(self.le_water_cut.text())
        self.surface_tension = float(self.le_surface_tension.text())
        self.N = float(self.le_N.text())
        self.Q = float(self.le_Q.text())
        self.stage_number = float(self.le_stage_number.text())
        self.pump = pump
        self.QBEM = QBEM
        self.catalogData = catalog_data
        self.df_sgle = None

        self.btn_calculate.clicked.connect(self.btn_calculate_clicked)
        self.btn_clear.clicked.connect(self.btn_clear_clicked)
        self.btn_exit.clicked.connect(self.btn_exit_clicked)
        self.btn_save.clicked.connect(self.btn_save_clicked)

        # plot initialize
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        plot_layout = QtWidgets.QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        self.widget_plot.setLayout(plot_layout)

    @QtCore.pyqtSlot()
    def btn_calculate_clicked(self):
        # update the input values
        self.oil_viscosity = float(self.le_oil_viscosity.text())
        self.oil_density = float(self.le_oil_density.text())
        self.water_viscosity = float(self.le_water_viscosity.text())
        self.water_density = float(self.le_water_density.text())
        self.water_cut = float(self.le_water_cut.text())
        self.surface_tension = float(self.le_surface_tension.text())
        self.N = float(self.le_N.text())
        self.Q = float(self.le_Q.text())
        self.stage_number = float(self.le_stage_number.text())
        self.VISE = 0.

        # return MainProgram.single_phase_calculation(self.pump,)
        if len(self.pump) < 5:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Critical)
            msg.setText("Please select the pump or input pump geometry!")
            msg.setWindowTitle("Input Error")
            msg.exec_()
        else:
            self.VISE = self.emulsion_rheology(self.oil_viscosity, self.water_viscosity, self.oil_density,
                                          self.water_density, self.water_cut, self.surface_tension,
                                          self.N, self.Q, self.stage_number)

            self.lbl_effective_viscosity.setText(str("%.5f" % self.VISE))

        # plot single-phase water data
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(self.catalogData['Q'], self.catalogData['P'], 'bo', markerfacecolor="None", label='Catalog')

        # plot calculated results
        QL = np.arange(500)[1:] * 50.0
        hpsgl = []
        hesgl = []
        heesgl = []
        hfsgl = []
        htsgl = []
        hresgl = []
        qlksgl = []
        self.pump['VISL'] = self.VISE
        self.pump['DENL'] = (1.-self.water_cut/100.)*self.oil_density + self.water_cut/100.*self.water_density
        for q in QL:
            HP_sgl, HE_sgl, HEE_sgl, HF_sgl, HT_sgl, HRE_sgl, QLK_sgl = MainProgram.single_phase_calculation(
                self.pump, self.QBEM, q)
            if HP_sgl > 0:
                hpsgl.append(HP_sgl)
                hesgl.append(HE_sgl)
                heesgl.append(HEE_sgl)
                hfsgl.append(HF_sgl)
                htsgl.append(HT_sgl)
                hresgl.append(HRE_sgl)
                qlksgl.append(QLK_sgl)
            else:
                break

        self.df_sgle = pd.DataFrame({'QL': QL[:len(hpsgl)], 'HP': hpsgl, 'HE': hesgl, 'HEE': heesgl, 'HF': hfsgl,
                                'HT': htsgl, 'HRE': hresgl, 'QLK': qlksgl})

        ax.plot(self.df_sgle.QL, self.df_sgle.HP, 'r-', label=r'$\mu_E$={}cp'.format(1000*round(self.VISE, 5)))
        ax.tick_params(direction='in')
        ax.set_title(r'Emulsion Performance, Water cut={}%'.format(self.water_cut))
        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel('P (psi)')
        ax.legend(frameon=False)
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def btn_clear_clicked(self):
        self.le_oil_viscosity.clear()
        self.le_oil_density.clear()
        self.le_water_viscosity.clear()
        self.le_water_density.clear()
        self.le_water_cut.clear()
        self.le_surface_tension.clear()
        self.le_N.clear()
        self.le_Q.clear()
        self.le_stage_number.clear()
        self.figure.clear()
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def btn_exit_clicked(self):
        msg = QtWidgets.QMessageBox.question(self, "Emulsion Module", "Do you want to quit emulsion module",
                                             QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel,
                                             QtWidgets.QMessageBox.Ok)
        if msg == QtWidgets.QMessageBox.Ok:
            self.accept()


    @QtCore.pyqtSlot()
    def btn_save_clicked(self):
        self.figure.savefig('Emulsion Figure.jpg', dpi=300)

        writer = pd.ExcelWriter('Output_Emulsion.xlsx', engine='xlsxwriter')
        self.df_sgle.to_excel(writer, sheet_name='Emulsion performance')
        writer.save()

        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The erosion calculation data and figure have been saved!")
        msg.setWindowTitle("Save Files")
        msg.exec_()

    # the model is based on Brinkman (1952) correlation and Zhang (2017, Fall) TUALP ABM
    def emulsion_rheology(self, VISO, VISW, DENO, DENW, WC, ST, N, Q, SN):
        """
        :param VISO:viscosity of oil (kg/m-s)
        :param VISW:viscosity of water (kg/m-s)
        :param DENO:density of oil (kg/m3)
        :param DENW:density of water (kg/m3)
        :param WC:  water cut (%)
        :param ST:  surface tension (N/m)
        :param N:   rotational speed (rpm)
        :param Q:   flow rate (bpd)
        :param SN:  stage number (-)
        :return:    effective viscosity (kg/m-s)
        """
        E = 3.   # exponential index
        Q = Q * 0.0000018401307361111
        f = N/60.
        WC = WC/100.
        miu_tilda = VISO/VISW
        phi_OI = miu_tilda**(1./E)/(1. + miu_tilda**(1./E))
        phi_WI = 1. - phi_OI
        phi_OE = 1. - (VISW/VISO)**(1./E)
        i = 0.

        # find the inversion point
        St = f*self.pump['VOI']/Q
        for i in np.arange(1000)/1000.:
            if i == 0:
                continue
            rouA = i * DENW + (1 - i) * DENO
            We = rouA * Q ** 2 / (ST * self.pump['VOI'])
            miu_M = VISW / (1 - (1-i) * phi_OE) ** E

            # assume oil in water
            Re = rouA*Q/(VISW*2*self.pump['R2'])
            C = (SN*We*Re)**0.15/(10*St**0.5)
            miu_E = VISW/(1 - (1-i))**E
            miu_A_OW = C*(miu_E - miu_M) + miu_M

            # assume water in oil
            Re = rouA * Q / (VISO * 2 * self.pump['R2'])
            C = (SN * We * Re) ** 0.15 / (10 * St ** 0.5)
            miu_E = VISO / (1 - i) ** E
            miu_A_WO = C * (miu_E - miu_M) + miu_M

            if np.abs(miu_A_OW-miu_A_WO)/miu_A_OW < 0.01:
                break

        if WC > i:
            # oil in water
            rouA = WC * DENW + (1 - WC) * DENO
            We = rouA * Q ** 2 / (ST * self.pump['VOI'])
            miu_M = VISW / (1 - (1 - WC) * phi_OE) ** E

            Re = rouA * Q / (VISW * 2 * self.pump['R2'])
            C = (SN * We * Re) ** 0.15 / (10 * St ** 0.5)
            miu_E = VISW / (1 - (1 - WC)) ** E
            miu_A = C * (miu_E - miu_M) + miu_M
        else:
            # water in oil
            rouA = WC * DENW + (1 - WC) * DENO
            We = rouA * Q ** 2 / (ST * self.pump['VOI'])
            miu_M = VISW / (1 - (1 - WC) * phi_OE) ** E

            Re = rouA * Q / (VISO * 2 * self.pump['R2'])
            C = (SN * We * Re) ** 0.15 / (10 * St ** 0.5)
            miu_E = VISO / (1 - WC) ** E
            miu_A = C * (miu_E - miu_M) + miu_M

        return miu_A


#######################################
# class for emulsion module calculation
#######################################
class DlgErosionModule(QtWidgets.QDialog, Ui_dlg_erosion_module):
    def __init__(self, pump, catalog_data, QBEM, parent=None):
        super(DlgErosionModule, self).__init__(parent)
        self.setupUi(self)

        # initialize new dialogue inputs
        self.le_RB.setValidator(QtGui.QDoubleValidator())
        self.le_RB.setText(str(0.004))
        self.le_NB.setValidator(QtGui.QDoubleValidator())
        self.le_NB.setText(str(6))
        self.le_LB.setValidator(QtGui.QDoubleValidator())
        self.le_LB.setText(str(0.01))
        self.le_SL1.setValidator(QtGui.QDoubleValidator())
        self.le_SL1.setText(str(0.000440))
        self.le_RLI1.setValidator(QtGui.QDoubleValidator())
        self.le_RLI1.setText(str(0.039403))
        self.le_RLO1.setValidator(QtGui.QDoubleValidator())
        self.le_RLO1.setText(str(0.01894))
        self.le_LG1.setValidator(QtGui.QDoubleValidator())
        self.le_LG1.setText(str(0.005))
        self.le_SL2.setValidator(QtGui.QDoubleValidator())
        self.le_SL2.setText(str(0.000330))
        self.le_RLI2.setValidator(QtGui.QDoubleValidator())
        self.le_RLI2.setText(str(0.04119))
        self.le_RLO2.setValidator(QtGui.QDoubleValidator())
        self.le_RLO2.setText(str(0.01894))
        self.le_LG2.setValidator(QtGui.QDoubleValidator())
        self.le_LG2.setText(str(0.005))
        self.le_HC1.setValidator(QtGui.QDoubleValidator())
        self.le_HC1.setText(str(0.002))
        self.le_HC2.setValidator(QtGui.QDoubleValidator())
        self.le_HC2.setText(str(0.002))
        self.le_HPW.setValidator(QtGui.QDoubleValidator())
        self.le_HPW.setText(str(0.001))
        self.le_QERO.setValidator(QtGui.QDoubleValidator())
        self.le_QERO.setText(str(3100))
        self.le_DENP.setValidator(QtGui.QDoubleValidator())
        self.le_DENP.setText(str(2600))
        self.le_DP.setValidator(QtGui.QDoubleValidator())
        self.le_DP.setText(str(0.0001))
        self.le_BH.setValidator(QtGui.QDoubleValidator())
        self.le_BH.setText(str(180))
        self.le_FS.setValidator(QtGui.QDoubleValidator())
        self.le_FS.setText(str(1))
        self.le_IMPW.setValidator(QtGui.QDoubleValidator())
        self.le_IMPW.setText(str(0.312))
        self.le_TT.setValidator(QtGui.QDoubleValidator())
        self.le_TT.setText(str(200))
        self.le_E.setValidator(QtGui.QDoubleValidator())
        self.le_E.setText(str(21E7))
        self.le_LST.setValidator(QtGui.QDoubleValidator())
        self.le_LST.setText(str(0.070358))
        self.le_DS.setValidator(QtGui.QDoubleValidator())
        self.le_DS.setText(str(0.0174752))
        self.le_TS1.setValidator(QtGui.QDoubleValidator())
        self.le_TS1.setText(str(0.00381))
        self.le_TS2.setValidator(QtGui.QDoubleValidator())
        self.le_TS2.setText(str(0.003175))
        self.le_DENPUMP.setValidator(QtGui.QDoubleValidator())
        self.le_DENPUMP.setText(str(7700))

        # read all inputs from dialogue
        self.RB = float(self.le_RB.text())
        self.NB = float(self.le_NB.text())
        self.LB = float(self.le_LB.text())
        self.SL1 = float(self.le_SL1.text())
        self.RLI1 = float(self.le_RLI1.text())
        self.RLO1 = float(self.le_RLO1.text())
        self.LG1 = float(self.le_LG1.text())
        self.SL2 = float(self.le_SL2.text())
        self.RLI2 = float(self.le_RLI2.text())
        self.RLO2 = float(self.le_RLO2.text())
        self.LG2 = float(self.le_LG2.text())
        self.HC1 = float(self.le_HC1.text())
        self.HC2 = float(self.le_HC2.text())
        self.HPW = float(self.le_HPW.text())
        self.QERO = float(self.le_QERO.text())
        self.DENP = float(self.le_DENP.text())
        self.DP = float(self.le_DP.text())
        self.BH = float(self.le_BH.text())
        self.FS = float(self.le_FS.text())
        self.IMPW = float(self.le_IMPW.text())
        self.TT = float(self.le_TT.text())
        self.E = float(self.le_E.text())
        self.LST = float(self.le_LST.text())
        self.DS = float(self.le_DS.text())
        self.TS1 = float(self.le_TS1.text())
        self.TS2 = float(self.le_TS2.text())
        self.DENPUMP = float(self.le_DENPUMP.text())
        self.EROM = self.comboBox_EROM.currentIndex()

        # read inputs from class initialization
        self.QBEM = QBEM * 0.15898 / 24.0 / 3600.0
        self.R1 = pump['R1']
        self.R2 = pump['R2']
        self.TB = pump['TB']
        self.TV = pump['TV']
        self.YI1 = pump['YI1']
        self.YI2 = pump['YI2']
        self.VOI = pump['VOI']
        self.VOD = pump['VOD']
        self.ASF = pump['ASF']
        self.ASB = pump['ASB']
        self.AB = pump['AB']
        self.AV = pump['AV']
        self.ADF = pump['ADF']
        self.ADB = pump['ADB']
        self.LI = pump['LI']
        self.LD = pump['LD']
        self.RLK = pump['RLK']
        self.LG = pump['LG']
        self.SL = pump['SL']
        self.EA = pump['EA']
        self.ZI = pump['ZI']
        self.ZD = pump['ZD']
        self.B1 = pump['B1']
        self.B2 = pump['B2']
        self.DENL = pump['DENL']
        self.VISL = pump['VISL']
        self.VISW = pump['VISW']
        self.N = pump['N']
        self.catalogData = catalog_data
        self.eroper = pd.DataFrame()
        self.allper = pd.DataFrame()

        # connect slot with signals
        self.btn_calculate.clicked.connect(self.btn_calculate_clicked)
        self.btn_clear.clicked.connect(self.btn_clear_clicked)
        self.btn_save.clicked.connect(self.btn_save_clicked)
        self.btn_exit.clicked.connect(self.btn_exit_clicked)
        self.radioButton_HQ.toggled.connect(self.radioButton_HQ_toggled)
        self.radioButton_HT.toggled.connect(self.radioButton_HT_toggled)
        self.radioButton_TS.toggled.connect(self.radioButton_TS_toggled)
        self.radioButton_LW.toggled.connect(self.radioButton_LW_toggled)
        self.tab2_stage.clicked.connect(self.tab2_stage_clicked)
        self.tab2_skirt.clicked.connect(self.tab2_skirt_clicked)
        self.tab2_balance.clicked.connect(self.tab2_balance_clicked)
        self.tab2_blade.clicked.connect(self.tab2_blade_clicked)

        # plot initialize
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.ax = self.figure.add_subplot(111)
        plot_layout = QtWidgets.QVBoxLayout()
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)
        self.tab1_Plot.setLayout(plot_layout)
        

    @QtCore.pyqtSlot()
    def btn_calculate_clicked(self):
        # initialize dataframe
        self.eroper = None
        self.allper = None

        # update inputs
        self.RB = float(self.le_RB.text())
        self.NB = float(self.le_NB.text())
        self.LB = float(self.le_LB.text())
        self.SL1 = float(self.le_SL1.text())
        self.RLI1 = float(self.le_RLI1.text())
        self.RLO1 = float(self.le_RLO1.text())
        self.LG1 = float(self.le_LG1.text())
        self.SL2 = float(self.le_SL2.text())
        self.RLI2 = float(self.le_RLI2.text())
        self.RLO2 = float(self.le_RLO2.text())
        self.LG2 = float(self.le_LG2.text())
        self.HC1 = float(self.le_HC1.text())
        self.HC2 = float(self.le_HC2.text())
        self.HPW = float(self.le_HPW.text())
        self.QERO = float(self.le_QERO.text()) * 0.15898 / 24.0 / 3600.0
        self.DENP = float(self.le_DENP.text())
        self.DP = float(self.le_DP.text())
        self.BH = float(self.le_BH.text())
        self.FS = float(self.le_FS.text())
        self.IMPW = float(self.le_IMPW.text())
        self.TT = float(self.le_TT.text())
        self.E = float(self.le_E.text())
        self.LST = float(self.le_LST.text())
        self.DS = float(self.le_DS.text())
        self.TS1 = float(self.le_TS1.text())
        self.TS2 = float(self.le_TS2.text())
        self.DENPUMP = float(self.le_DENPUMP.text())
        self.EROM = self.comboBox_EROM.currentIndex() + 1       # convert comboBox index into Fortran index

        if self.TT > 1e4:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Warning)
            msg.setText("The time input exceeds maximum value!")
            msg.setWindowTitle("Warning")
            msg.exec_()
        elif self.DP > self.SL1 and self.DP > self.SL2:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Warning)
            msg.setText("The sand particle size is too big!")
            msg.setWindowTitle("Warning")
            msg.exec_()
        else:
            # call subroutine to calculate
            _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, EROPER, ALLPER = \
            ESPERO.ero_mod.eroper_mod(self.QBEM, self.R1,   self.R2,   self.N,   self.B1, self.B2, self.TB,   self.TV,
                                      self.YI1,  self.YI2,  self.VOI,  self.VOD, self.AB, self.AV, self.ASF,  self.ASB,
                                      self.ADF,  self.ADB,  self.EA,   self.LI,  self.LD, self.ZI, self.ZD,   self.DENL,
                                      self.VISL, self.VISW, self.RB,   self.NB,  self.LB, self.SL1,self.RLI1, self.RLO1,
                                      self.LG1,  self.SL2,  self.RLI2, self.RLO2,self.LG2,self.HC1,self.HC2,  self.HPW,
                                      self.EROM, self.QERO, self.DENP, self.DP,  self.BH, self.FS, self.IMPW, self.TT,
                                      self.E,    self.LST,  self.DS,   self.TS1, self.TS2,self.DENPUMP)
            self.allper = pd.DataFrame(ALLPER)
            self.allper = self.allper.loc[:, (self.allper != 0).any(axis=0)]        # drop zero columns
            self.allper = self.allper[~(self.allper == 0).any(axis=1)]              # drop zero rows
            self.allper = self.allper[self.allper.columns[-2:]]                     # get last two columns
            self.allper.columns = ['Q', 'P']
            self.eroper = pd.DataFrame({'t': EROPER[:, 0], 'HP': EROPER[:, 1], 'SL12': EROPER[:, 2], 'SL22': EROPER[:, 3],
                                        'TS11': EROPER[:, 4], 'TS12': EROPER[:, 5], 'TBR': EROPER[:, 6], 'Life': EROPER[:, 7], 'OKA': EROPER[:, 8],
                                        'DNV': EROPER[:, 9], 'AHLERT': EROPER[:, 10], 'MANSOURI': EROPER[:, 11], 'ZHANG': EROPER[:, 12],
                                        'HAUGEN': EROPER[:, 13]})
            # print self.TT

            # plot single-phase water data
            self.tabWidget.setCurrentIndex(0)   # set current tab bar
            self.figure.clear()
            df = self.eroper[['t', 'HP']]
            df = df[~(df == 0).any(axis=1)]
            ax = self.figure.add_subplot(111)
            ax.plot(df['t'], df['HP'], 'r-')
            if df.iloc[-1]['t'] < self.TT:
                ax.annotate(r'Pump fails at %d day' % (df.iloc[-1]['t']), xy=(df.iloc[-1]['t'], df.iloc[-1]['HP']), xytext=(df.iloc[-1]['t']/2, (df.iloc[1]['HP']+df.iloc[-1]['HP'])/2),
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1))
            else:
                ax.annotate(r'Pump survives!', xy=(df.iloc[-1]['t'], df.iloc[-1]['HP']), xytext=(df.iloc[-1]['t']/1.5, (df.iloc[1]['HP']+df.iloc[-1]['HP'])/2),
                    arrowprops=dict(facecolor='black', shrink=0.1, width=1))
            ax.tick_params(direction='in')
            ax.set_title(r'Pump head at production rate versus erosion time')
            ax.set_xlabel(r'Time (day)')
            ax.set_ylabel('P (psi)')
            self.canvas.draw()

    @QtCore.pyqtSlot()
    def radioButton_HQ_toggled(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(self.catalogData['Q'], self.catalogData['P'], 'bo-', markerfacecolor="None", label='Catalog')
        ax.plot(self.allper.Q, self.allper.P, 'r-', label="Pump performance after erosion")

        ax.tick_params(direction='in')
        ax.set_title(r'Sand weight concentration={}%'.format(self.HPW))
        ax.set_xlabel(r'$Q_L$ (bpd)')
        ax.set_ylabel('P (psi)')
        ax.legend(frameon=False)
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def radioButton_HT_toggled(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.clear()
        df = self.eroper[['t', 'HP']]
        df = df[~(df == 0).any(axis=1)]
        ax = self.figure.add_subplot(111)
        ax.plot(df['t'], df['HP'], 'r-')
        if df.iloc[-1]['t'] < self.TT:
            ax.annotate(r'Pump fails at %d day' % (df.iloc[-1]['t']), xy=(df.iloc[-1]['t'], df.iloc[-1]['HP']), xytext=(df.iloc[-1]['t']/2, (df.iloc[1]['HP']+df.iloc[-1]['HP'])/2),
                arrowprops=dict(facecolor='black', shrink=0.1, width=1))
        else:
            ax.annotate(r'Pump survives!', xy=(df.iloc[-1]['t'], df.iloc[-1]['HP']), xytext=(df.iloc[-1]['t']/1.5, (df.iloc[1]['HP']+df.iloc[-1]['HP'])/2),
                arrowprops=dict(facecolor='black', shrink=0.1, width=1))
        ax.tick_params(direction='in')
        ax.set_title(r'Pump head at production rate versus erosion time')
        ax.set_xlabel(r'Time (day)')
        ax.set_ylabel('P (psi)')
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def radioButton_TS_toggled(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        df = self.eroper[['Life', 'OKA', 'DNV', 'AHLERT', 'MANSOURI', 'ZHANG', 'HAUGEN']]
        df = df[~(df == 0).any(axis=1)]
        ax.plot(df['Life'], df['OKA'], color='red', label='OKA')
        ax.plot(df['Life'], df['DNV'], color='blue', label='DNV')
        ax.plot(df['Life'], df['AHLERT'], color='black', label='AHLERT')
        ax.plot(df['Life'], df['MANSOURI'], color='green', label='MANSOURI')
        ax.plot(df['Life'], df['ZHANG'], color='purple', label='ZHANG')
        ax.plot(df['Life'], df['HAUGEN'], color='magenta', label='HAUGEN')
        ax.tick_params(direction='in')
        ax.legend(frameon=False, loc=10, fontsize=7)
        ax.set_title(r'Pump Failure Possibility')
        ax.set_xlabel(r'Erosion time (day)')
        ax.set_ylabel(r'Possibility (%)')
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def radioButton_LW_toggled(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        df = self.eroper[['t', 'SL12', 'SL22', 'TS12', 'TS11', 'TBR']]
        df = df[~(df == 0).any(axis=1)]
        ax.plot(df['t'], df['SL12'], color='red', label='Skirt gap width')
        ax.plot(df['t'], df['SL22'], color='blue', label='Balance gap width')
        ax.plot(df['t'], df['TS11'], color='black', label='Skirt ring thickness')
        ax.plot(df['t'], df['TS12'], color='green', label='Balance ring thickness')
        ax.plot(df['t'], df['TBR'], color='purple', label='Impeller blade thickness')
        ax.tick_params(direction='in')
        ax.legend(frameon=False, fontsize=7)
        ax.set_title(r'Level of Wear')
        ax.set_ylabel(r'Width/thickness (mm)')
        ax.set_xlabel(r'Erosion time (day)')
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def btn_clear_clicked(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.clear()
        self.canvas.draw()

    @QtCore.pyqtSlot()
    def btn_save_clicked(self):
        self.tabWidget.setCurrentIndex(0)  # set current tab bar
        self.figure.savefig('Erosion Figure.jpg', dpi=300)

        writer = pd.ExcelWriter('Output_Erosion.xlsx', engine='xlsxwriter')
        self.allper.to_excel(writer, sheet_name='Erosion performance curve')
        self.eroper.to_excel(writer, sheet_name='Erosion analysis')
        writer.save()

        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setText("The erosion calculation data and figure have been saved!")
        msg.setWindowTitle("Save Files")
        msg.exec_()

    @QtCore.pyqtSlot()
    def btn_exit_clicked(self):
        msg = QtWidgets.QMessageBox.question(self, "Erosion Module", "Do you want to quit erosion module",
                                             QtWidgets.QMessageBox.Ok | QtWidgets.QMessageBox.Cancel,
                                             QtWidgets.QMessageBox.Ok)
        if msg == QtWidgets.QMessageBox.Ok:
            self.accept()

    @QtCore.pyqtSlot()
    def tab2_stage_clicked(self):
        self.tab2_frame.setStyleSheet("background-color: rgb(255, 255, 255);\n"
                                      "border-image: url(:/esp_pics/Stage.png);")

    @QtCore.pyqtSlot()
    def tab2_skirt_clicked(self):
        self.tab2_frame.setStyleSheet("background-color: rgb(255, 255, 255);\n"
                                      "border-image: url(:/esp_pics/Skirt.png);")

    @QtCore.pyqtSlot()
    def tab2_balance_clicked(self):
        self.tab2_frame.setStyleSheet("background-color: rgb(255, 255, 255);\n"
                                      "border-image: url(:/esp_pics/Balance.png);")

    @QtCore.pyqtSlot()
    def tab2_blade_clicked(self):
        self.tab2_frame.setStyleSheet("background-color: rgb(255, 255, 255);\n"
                                      "border-image: url(:/esp_pics/Blade.png);")

######################
# run the main program
######################
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    ui = MainProgram()
    ui.show()
    sys.exit(app.exec_())
