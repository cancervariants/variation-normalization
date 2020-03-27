import { Component } from "react";
import { ValidationResponse, ValidationResult } from "../types/VarlexTypes";
import { Divider, Input, InputOnChangeData, Table, Segment, Button, Icon, Select } from "semantic-ui-react";
import React from "react";
import { getValidations } from "../services/VarlexApi";
import { ValidationResultRow, CheckboxCallback } from "./ValidationResultRow";

type State = {
    validationResponse: ValidationResponse | null;
    activeTimeout: number | null;
    selectedResults: Set<ValidationResult>;
}

export class ValidationTable extends Component<{}, State> {
    state: State = {
        validationResponse: null,
        activeTimeout: null,
        selectedResults: new Set()
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
        let keyCount = 0
        if (this.state.validationResponse) {
            const summary = this.state.validationResponse.validationSummary
            components.push(<Divider key={keyCount++} />)

            if (summary.validResults.length > 0 || summary.invalidResults.length > 0) {
                components.push(this.headerSegment(keyCount++))
            }

            if (summary.validResults.length > 0) {
                components.push(this.resultTable(summary.validResults, keyCount++))
                components.push(<Divider key={keyCount++} />)
            }
            if (summary.invalidResults.length > 0) {
                components.push(this.resultTable(summary.invalidResults, keyCount++))
                components.push(<Divider key={keyCount++} />)
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

    private headerSegment = (index: number): JSX.Element => {

        if (this.state.selectedResults.size === 0) {
            return <Segment secondary key={index}>Please select variants for export</Segment>
        } else {
            let formatOpts = [
                { key: 'json', text: 'JSON', value: 'json' },
                { key: 'vrs', text: 'VRS', value: 'vrs' },
            ]

            return (
                <Segment.Group key={index}>
                    <Segment secondary>
                        <p>{this.state.selectedResults.size} variant(s) selected.</p>
                    </Segment>
                    <Segment clearing secondary>
                        <Select options={formatOpts} defaultValue='vrs' />
                        <Button animated floated="right">
                            <Button.Content visible>Export</Button.Content>
                            <Button.Content hidden>
                                <Icon name='external share' />
                            </Button.Content>
                        </Button>
                    </Segment>
                </Segment.Group>
            )
        }
    }

    private resultTable = (content: ValidationResult[], index: number): JSX.Element => {
        const color = content.length > 0 && content[0].isValid ? 'green' : 'red'
        return (
            <Table color={color} key={index}>
                {this.tableHeader(content[0].isValid)}
                {this.tableContents(content)}
            </Table>
        );
    }

    private tableHeader = (valid: boolean): JSX.Element => {
        if (valid) {
            return (<Table.Header>
                <Table.Row>
                    <Table.HeaderCell></Table.HeaderCell>
                    <Table.HeaderCell>Identified Variant</Table.HeaderCell>
                    <Table.HeaderCell>Reference Build</Table.HeaderCell>
                    <Table.HeaderCell>Classification</Table.HeaderCell>
                    <Table.HeaderCell>Description</Table.HeaderCell>
                </Table.Row>
            </Table.Header>);
        } else {
            return (<Table.Header>
                <Table.Row>
                    <Table.HeaderCell></Table.HeaderCell>
                    <Table.HeaderCell>Identified Variant</Table.HeaderCell>
                    <Table.HeaderCell>Reference Build</Table.HeaderCell>
                    <Table.HeaderCell>Classification</Table.HeaderCell>
                    <Table.HeaderCell>Description</Table.HeaderCell>
                    <Table.HeaderCell>Validation Errors</Table.HeaderCell>
                </Table.Row>
            </Table.Header>);
        }
    }

    private tableContents = (content: ValidationResult[]): JSX.Element => {

        const rows = content.map((r: ValidationResult, index: number) =>
            <ValidationResultRow result={r} key={index} index={index} onSelected={this.onResultSelected} />
        );
        return <Table.Body>{rows}</Table.Body>;
    }


    private onResultSelected: CheckboxCallback = (res: ValidationResult, newVal: boolean): void => {
        this.setState((prevState) => {
            return {
                ...prevState,
                selectedResult: newVal ? prevState.selectedResults.add(res) : prevState.selectedResults.delete(res)
            }
        })
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