import React from 'react';
import './App.css';
import { TokenTable } from './components/TokenTable';
import { Menu, Container, Tab } from 'semantic-ui-react';
import { ClassificationTable } from './components/ClassificationTable';
import { SwaggerDocs } from './components/SwaggerDocs';
import { ValidationTable } from './components/ValidationTable';

const App: React.FC = () => {
  const tabPanes = [
    { menuItem: 'Tokenizer', render: () => <TokenTable /> },
    { menuItem: 'Classifier', render: () => <ClassificationTable /> },
    { menuItem: 'Validator', render: () => <ValidationTable /> },
    { menuItem: 'OpenAPI Docs', render: () => <SwaggerDocs /> },
  ]
  return (
    <div>
      <Menu fixed='top' inverted>
        <Container>
          <Menu.Item header>Varlex Prototype</Menu.Item>
          <Menu.Item as='a' href="/openapi.json">OpenAPI JSON</Menu.Item>
        </Container>
      </Menu>

      <Container style={{ marginTop: '7em' }}>
        <Tab panes={tabPanes} menu={{ pointing: true }} />
      </Container>
    </div >
  );
}

export default App;
