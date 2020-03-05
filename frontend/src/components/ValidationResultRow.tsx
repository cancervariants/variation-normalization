import React from "react";
import { Component } from "react";
import { ValidationResult } from "../types/VarlexTypes";
import { Table, Checkbox } from "semantic-ui-react";


export type CheckboxCallback = (res: ValidationResult, newVal: boolean) => void

type State = {
    selected: boolean
}

type Props = {
    result: ValidationResult;
    index: number;
    onSelected: CheckboxCallback;
}

export class ValidationResultRow extends Component<Props, State> {

    state: State = {
        selected: false
    }

    render() {
        return this.tableRow(this.props.result, this.props.index)
    }


    private onCheckboxClicked = (): void => {
        this.props.onSelected(this.props.result, !this.state.selected);
        this.setState((prevState) => {
            return {
                ...prevState,
                selected: !prevState.selected
            }
        });
    }


    private tableRow = (res: ValidationResult, index: number): JSX.Element => {
        if (res.isValid) {
            return (<Table.Row key={index} active={this.state.selected}>
                <Table.Cell><Checkbox checked={this.state.selected} onChange={this.onCheckboxClicked} /></Table.Cell>
                <Table.Cell>{res.conciseDescription}</Table.Cell>
                <Table.Cell>hg38</Table.Cell>
                <Table.Cell>{res.classification.classificationType}</Table.Cell>
                <Table.Cell>{res.humanDescription}</Table.Cell>
            </Table.Row>)
        } else {
            return (<Table.Row key={index} active={this.state.selected}>
                <Table.Cell><Checkbox checked={this.state.selected} onChange={this.onCheckboxClicked} /></Table.Cell>
                <Table.Cell>{res.conciseDescription}</Table.Cell>
                <Table.Cell>hg38</Table.Cell>
                <Table.Cell>{res.classification.classificationType}</Table.Cell>
                <Table.Cell>{res.humanDescription}</Table.Cell>
                <Table.Cell>{res.errors.join(', ')}</Table.Cell>
            </Table.Row>)
        }
    }
}