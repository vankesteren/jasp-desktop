//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//
import QtQuick          2.8
import QtQuick.Layouts  1.3
import JASP.Controls    1.0
import JASP.Theme		1.0
import JASP.Widgets     1.0

Form
{

    VariablesForm
    {
        height: 420
        AvailableVariablesList
        {
            name: "availableVariables"
        }
        AssignedVariablesList
        {
            name: "variables"
            title: qsTr("Variables")
            suggestedColumns: ["scale"]
            id: variables
        }
        AssignedVariablesList
        {
            name: "regressions"
            title: qsTr("Regressions")
            suggestedColumns: ["scale"]
            allowedColumns: ["scale"]
        }
        AssignedVariablesList
        {
            name: "categorical"
            id: categorical
            title: qsTr("Factor")
            suggestedColumns: ["ordinal", "nominal"]
            allowedColumns: ["ordinal", "nominal", "nominalText"]
        }
        AssignedVariablesList
        {
            name: "covariates"
            title: qsTr("Time-varying covariates")
            suggestedColumns: ["scale"]
            debug: true
            //listViewType: "RepeatedMeasures"
            //source: variables.name
        }
    }

    AssignedVariablesList
    {
        name: "timings"
        title: qsTr("Timings")
        source: "variables"
        draggable: false
        preferredWidth: form.width * 2 / 5
        preferredHeight: 140

        rowComponents:
        [
               Component { DoubleField {
                    name: "timing"
                    negativeValues: true
                    defaultValue: rowIndex
                    useExternalBorder: false
                } }
        ]
    }

    Section
    {
        title: qsTr("Model Options")
        Group
        {
            title: qsTr("Growth curve shape")
            CheckBox
            {
                label: qsTr("Intercept")
                name: "intercept"
                id: intercept
                checked: true
            }
            CheckBox
            {
                label: qsTr("Linear")
                name: "linear"
                id: linear
                enabled: intercept.checked
                checked: true
            }
            CheckBox
            {
                label: qsTr("Quadratic")
                name: "quadratic"
                id: quadratic
                enabled: linear.checked
                checked: if (!linear.checked) return false
            }
            CheckBox
            {
                label: qsTr("Cubic")
                name: "cubic"
                id: cubic
                enabled: quadratic.checked
                checked: if (!quadratic.checked) return false
            }
        }

        Group
        {
            CheckBox
            {
                label: qsTr("Covarying latent curve")
                name: "covar"
                enabled: linear.checked
                checked: true
            }
        }
    }

    Section
    {
        title: qsTr("Additional Output")
        GroupBox
        {
            CheckBox { label: qsTr("Additional Fit Measures")   ; name: "outputAdditionalFitMeasures"   }
            CheckBox { label: qsTr("R-Squared")                 ; name: "rsquared"                      }
            CheckBox { label: qsTr("Standardized estimates")    ; name: "std"                           }
        }
        GroupBox
        {
            CheckBox { label: qsTr("Implied covariance matrix")  ; name: "impliedCov" }
            CheckBox { label: qsTr("Residual covariance matrix") ; name: "residCov"   }
            CheckBox { label: qsTr("Show lavaan syntax")         ; name: "showSyntax" }
        }
    }

    Section
    {
        text: qsTr("Plots")
        GroupBox
        {
            title: "Plots"
            DropDown
            {
                name: "colorPalette"
                label: qsTr("Color palette")
                indexDefaultValue: 0
                values:
                [
                    { label: qsTr("Colorblind"),		value: "colorblind"		},
                    { label: qsTr("Colorblind Alt."),	value: "colorblind3"	},
                    { label: qsTr("Viridis"),			value: "viridis"		},
                    { label: qsTr("ggplot2"),			value: "ggplot2"		},
                    { label: qsTr("Gray"),				value: "gray"			}
                ]
            }
            CheckBox {
                text: qsTr("Curve plot")
                name: "curveplot"
                DropDown {
                    name: "plot_categorical"
                    label: "Colour lines by"
                    source: "categorical"
                    addEmptyValue: true
                    enabled: categorical.count > 0
                }

                IntegerField {
                    name: "plot_max_n"
                    text: "Maximum number of lines"
                    defaultValue: 150
                    negativeValues: false
                }
            }
            CheckBox { text: qsTr("Misfit plot")    ; name: "misfitplot"   ; debug: true }
            CheckBox
            {
                text: qsTr("Model plot")
                name: "pathplot"
                CheckBox { text: qsTr("Show parameters") ; name: "plotpars"  }
                CheckBox { text: qsTr("Show means")      ; name: "plotmeans" }
            }
        }
    }

    Section
    {
        text: qsTr("Advanced")
        RadioButtonGroup
        {
            title: qsTr("Emulation")
            name: "mimic"
            RadioButton { text: qsTr("None")  ; name: "lavaan"  ; checked: true }
            RadioButton { text: qsTr("Mplus") ; name: "Mplus" }
            RadioButton { text: qsTr("EQS")   ; name: "EQS"   }
        }

        RadioButtonGroup {
            title: qsTr("Missing value handling")
            name: "missing"
            RadioButton { text: qsTr("Full Information Maximum Likelihood") ; name: "fiml" ; checked: true }
            RadioButton { text: qsTr("Exclude cases listwise")              ; name: "listwise"             }
        }

        GroupBox
        {
            title: qsTr("Error calculation")
            CIField { text: qsTr("CI width"); name: "ciWidth" }
            RadioButtonGroup
            {
                title: qsTr("Method")
                name: "se"
                RadioButton { text: qsTr("Standard")  ; name: "standard" ; checked: true }
                RadioButton { text: qsTr("Robust")    ; name: "robust" }
                RadioButton {
                    text: qsTr("Bootstrap CI")
                    name: "bootstrap"
                    IntegerField {
                        text: qsTr("Bootstrap samples")
                        name: "bootstrapNumber"
                        defaultValue: 1000
                        min: 100
                        max: 1000000
                    }
                }
            }
        }

        RadioButtonGroup
        {
            title: qsTr("Estimator")
            name: "estimator"
            RadioButton { text: qsTr("Auto") ; name: "default"; checked: true }
            RadioButton { text: qsTr("ML")   ; name: "ML"       }
            RadioButton { text: qsTr("GLS")  ; name: "GLS"      }
            RadioButton { text: qsTr("WLS")  ; name: "WLS"      }
            RadioButton { text: qsTr("ULS")  ; name: "ULS"      }
            RadioButton { text: qsTr("DWLS") ; name: "DWLS"     }
        }

        GroupBox
        {
            title: qsTr("Options")
            debug: true
            CheckBox { text: qsTr("Fix manifest intercepts to zero") ; name: "fixManifestInterceptsToZero" }
            CheckBox { text: qsTr("Fix latent intercepts to zero")   ; name: "fixLatentInterceptsToZero"   ; checked: true }
            CheckBox { text: qsTr("Omit residual single indicator")  ; name: "omitResidualSingleIndicator" ; checked: true }
            CheckBox { text: qsTr("Residual variances")              ; name: "residualVariances"           ; checked: true }
            CheckBox { text: qsTr("Correlate exogenous latents")     ; name: "correlateExogenousLatents"   ; checked: true }
            CheckBox { text: qsTr("Add thresholdds")                 ; name: "addThresholds"               ; checked: true }
            CheckBox { text: qsTr("Add scalings parameters")         ; name: "addScalingParameters"        ; checked: true }
            CheckBox { text: qsTr("Correlate dependent variables")   ; name: "correlateDependentVariables" ; checked: true }
        }
    }    

    Section
    {
        text: qsTr("Multigroup LGCM")
        debug: true
        DropDown
        {
            label: qsTr("Grouping variable") ;
            name: "groupvar";
            showVariableTypeIcon: true;
            addEmptyValue: true;
        } // No model or syncModels: it takes all variables per default
    }
}


