import { Component } from "react";
import { ValidationResponse, ValidationResult } from "../types/VarlexTypes";
import { Divider, Input, InputOnChangeData, Table } from "semantic-ui-react";
import React from "react";
import { getValidations } from "../services/VarlexApi";

type State = {
    validationResponse: ValidationResponse | null;
    activeTimeout: number | null;

}

export class ValidationTable extends Component<{}, State> {
    state: State = {
        validationResponse: null,
        activeTimeout: null
    }

    render() {
        return (
            <div>
                <Divider />
                <h3>Validation Testing</h3>
                <Input icon='search' placeholder='Validate' onChange={this.onSearchChanged} />
                {this.mainBody()}
            </div >
        );
    }

    private mainBody = (): JSX.Element => {
        let components = [];
        if (this.state.validationResponse) {
            const summary = this.state.validationResponse.validationSummary
            if (summary.validResults.length > 0) {
                components.push(<Divider />)
                components.push(this.resultTable(summary.validResults))
            }
            if (summary.invalidResults.length > 0) {
                components.push(<Divider />)
                components.push(this.resultTable(summary.invalidResults))
            }
        }

        if (components.length === 0) {
            return <div>No Matches...</div>;
        } else {
            return <div>
                {components}
            </div>
        }
    }

    private resultTable = (content: ValidationResult[]): JSX.Element => {
        const color = content.length > 0 && content[0].isValid ? 'green' : 'red'
        return (
            <Table color={color} key={color}>
                {this.tableHeader()}
                {this.tableContents(content)}
            </Table>
        );
    }

    private tableHeader = (): JSX.Element => {
        return (<Table.Header>
            <Table.Row>
                <Table.HeaderCell>Identified Variant</Table.HeaderCell>
                <Table.HeaderCell>Classification</Table.HeaderCell>
                <Table.HeaderCell>Description</Table.HeaderCell>
                <Table.HeaderCell>Validation Errors</Table.HeaderCell>
            </Table.Row>
        </Table.Header>);
    }

    private tableContents = (content: ValidationResult[]): JSX.Element => {
        const rows = content.map((r: ValidationResult, index: number) =>
            <Table.Row key={index}>
                <Table.Cell>{r.conciseDescription}</Table.Cell>
                <Table.Cell>{r.classification.classificationType}</Table.Cell>
                <Table.Cell>{r.humanDescription}</Table.Cell>
                <Table.Cell>{r.errors.join(', ')}</Table.Cell>
            </Table.Row>
        );
        return <Table.Body>{rows}</Table.Body>;
    }

    private onSearchChanged = (event: React.ChangeEvent<HTMLInputElement>, data: InputOnChangeData) => {
        if (this.state.activeTimeout) {
            window.clearTimeout(this.state.activeTimeout);
        }

        let searchTerm = event.target.value;
        let newTimer = window.setTimeout(() => { this.validate(searchTerm) }, 500);

        this.setState((prevState) => {
            return {
                ...prevState,
                activeTimeout: newTimer
            }
        });
    }

    private validate = (searchTerm?: string) => {
        getValidations(searchTerm || '')
            .then(validationResponse => this.setState({ validationResponse: validationResponse }))
    }
}