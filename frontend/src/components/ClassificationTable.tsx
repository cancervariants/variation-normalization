import React, { Component } from "react";
import { Table, Input, InputOnChangeData, Divider } from "semantic-ui-react";
import { ClassificationResponse, Classification } from "../types/VarlexTypes";
import { getClassifications } from "../services/VarlexApi"

type State = {
    classificationResponse: ClassificationResponse | null;
    activeTimeout: number | null;
};

export class ClassificationTable extends Component<{}, State> {

    state: State = {
        classificationResponse: null,
        activeTimeout: null
    }

    render() {
        return (
            <div>
                <Divider />
                <h3>Classification Testing</h3>
                <Input icon='search' placeholder='Classify' onChange={this.onSearchChanged} />
                <Table celled>
                    <Table.Header>
                        <Table.Row>
                            <Table.HeaderCell>Classification</Table.HeaderCell>
                            <Table.HeaderCell>Confidence</Table.HeaderCell>
                            <Table.HeaderCell>Matched Tokens</Table.HeaderCell>
                            <Table.HeaderCell>Unmatched Tokens</Table.HeaderCell>
                        </Table.Row>
                    </Table.Header>
                    {this.tableContents()}
                </Table>
            </div>
        );
    }

    private onSearchChanged = (event: React.ChangeEvent<HTMLInputElement>, data: InputOnChangeData) => {
        if (this.state.activeTimeout) {
            window.clearTimeout(this.state.activeTimeout);
        }

        let searchTerm = event.target.value;
        let newTimer = window.setTimeout(() => { this.classify(searchTerm) }, 500);

        this.setState((prevState) => {
            return {
                ...prevState,
                activeTimeout: newTimer
            }
        });
    }

    private classify = (searchTerm?: string) => {
        getClassifications(searchTerm || '')
            .then(classificationResponse => this.setState({ classificationResponse: classificationResponse }));
    }

    private tableContents = (): JSX.Element => {
        if (this.state.classificationResponse && this.state.classificationResponse.classifications.length > 0) {
            const rows = this.state.classificationResponse.classifications.map((c: Classification, index: number) =>
                <Table.Row key={index}>
                    <Table.Cell>{c.classificationType}</Table.Cell>
                    <Table.Cell>{c.confidence}</Table.Cell>
                    <Table.Cell>{c.matchingTokens.join(', ')}</Table.Cell>
                    <Table.Cell>{c.nonMatchingTokens.join(', ')}</Table.Cell>
                </Table.Row>
            );
            return <Table.Body>{rows}</Table.Body>
        } else {
            return <Table.Body><Table.Row><Table.Cell>No Matches...</Table.Cell></Table.Row></Table.Body>;
        }
    }
}
