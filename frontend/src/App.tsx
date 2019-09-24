import React from 'react';
import './App.css';
import { TokenTable } from './components/TokenTable';
import { Menu, Container } from 'semantic-ui-react';

const App: React.FC = () => {
  return (
    <div>
      <Menu fixed='top' inverted>
        <Container>
          <Menu.Item header>
            Varlex Prototype
        </Menu.Item>
          <Menu.Item as='a'>
            Tokenization Testing
        </Menu.Item>
        </Container>
      </Menu>

      <Container style={{ marginTop: '7em' }}>
        <TokenTable />
      </Container>
    </div >
  );
}

export default App;
