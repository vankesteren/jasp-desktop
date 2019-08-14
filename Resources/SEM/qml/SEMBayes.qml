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
import QtQuick 2.8
import QtQuick.Layouts 1.3
import JASP.Controls 1.0

Form
{
    usesJaspResults: true
	columns: 1
	TextArea
	{
		title: qsTr("Enter lavaan syntax below")
		name: "model"
		textType: "lavaan"
	}
    Section
    {
        title: qsTr("Estimation")
        GroupBox {
            IntegerField
            {
                name: "burnin"
                label: qsTr("Burn-in iterations")
                defaultValue: 500
                min: 500
            }
            IntegerField
            {
                name: "nsamples"
                label: qsTr("Samples per chain")
                defaultValue: 1000
                min: 500
            }
            IntegerField
            {
                name: "nchains"
                label: qsTr("Number of chains")
                defaultValue: 2
                min: 1
                max: 12
            }
        }
    }
    Section
    {
        title: qsTr("Diagnostics")
        GroupBox
        {
            CheckBox { label: qsTr("Trace plots")           ; name: "traceplot" }
            CheckBox { label: qsTr("Parameter diagnostics") ; name: "pardiag"  }
        }
    }
}
