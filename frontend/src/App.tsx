import React from 'react';
import './App.css';
import { TokenTable } from './components/TokenTable';
import { Menu, Container, Divider } from 'semantic-ui-react';
import { ClassificationTable } from './components/ClassificationTable';
import { SwaggerDocs } from './components/SwaggerDocs';

const App: React.FC = () => {
  return (
    <div>
      <Menu fixed='top' inverted>
        <Container>
          <Menu.Item header>Varlex Prototype</Menu.Item>
          <Menu.Item as='a' href="/openapi.json">OpenAPI JSON</Menu.Item>
        </Container>
      </Menu>

      <Container style={{ marginTop: '7em' }}>
        <TokenTable />
        <Divider section />
        <ClassificationTable />
        <Divider section />
        <SwaggerDocs />
      </Container>
    </div >
  );
}

export default App;
