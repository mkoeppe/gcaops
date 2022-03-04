r"""
Graph file view
"""
from abc import ABC, abstractmethod
import sqlite3
from math import ceil
from gcaops.util.undirected_graph_sage import undirected_graph_from_encoding, undirected_graph_to_encoding
from gcaops.util.directed_graph_sage import directed_graph_from_encoding, directed_graph_to_encoding
from gcaops.util.formality_graph_sage import formality_graph_from_encoding, formality_graph_to_encoding

class GraphFileView(ABC):
    """
    Graph database file view
    """
    def __init__(self, filename, num_vertices, num_edges):
        """
        Initialize this graph database file view.

        INPUT:

        - ``filename`` -- a string, the path to an SQLite database file (the file will be created if it does not yet exist)

        - ``num_vertices`` -- a natural number, the number of vertices in each graph

        - ``num_edges`` -- a natural number, the number of edges in each graph
        """
        self._filename = filename
        self._num_vertices = num_vertices
        self._num_edges = num_edges
        self._con = sqlite3.connect(filename)
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(name) FROM sqlite_master WHERE type = "table" AND name = "graphs"')
        if result.fetchone()[0] == 0:
            cur.execute('CREATE TABLE graphs (id INTEGER PRIMARY KEY, graph VARCHAR({0:d}))'.format(self._encoding_length()))
            self._con.commit()

    @abstractmethod
    def _encoding_to_graph(self, enc):
        """
        Return the graph corresponding to the given plain text encoding (e.g. in the graph6 or digraph6 format).
        """
        pass

    @abstractmethod
    def _graph_to_encoding(self, g):
        """
        Return the plain text encoding of the given graph (e.g. in the graph6 or digraph6 format).
        """
        pass

    @abstractmethod
    def _encoding_length(self):
        """
        Return the length of the plain text encoding of each graph in this database file.
        """
        pass

    def __getitem__(self, index):
        """
        Return the graph at the given index in this database file.
        """
        cur = self._con.cursor()
        result = cur.execute('SELECT graph FROM graphs WHERE id = ?', (int(index) + 1,)).fetchone()
        if result is not None:
            return self._encoding_to_graph(result[0])
        else:
            raise IndexError

    def index(self, g):
        """
        Return the index of the given graph in this database file.
        """
        encoding = self._graph_to_encoding(g)
        cur = self._con.cursor()
        result = cur.execute('SELECT id FROM graphs WHERE graph = ?', (encoding,)).fetchone()
        if result is not None:
            return result[0] - 1
        else:
            raise ValueError

    def __iter__(self):
        """
        Return an iterator over the graphs in this database file.
        """
        cur = self._con.cursor()
        for row in cur.execute('SELECT graph FROM graphs'):
            yield self._encoding_to_graph(row[0])

    def __len__(self):
        """
        Return the number of graphs in this database file.
        """
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(*) FROM graphs').fetchone()
        return result[0]

    def append(self, g):
        """
        Insert the given graph into this database file.
        """
        cur = self._con.cursor()
        cur.execute('INSERT INTO graphs (graph) VALUES (?)', (self._graph_to_encoding(g),))

    def commit(self):
        """
        Commit any changes to this database file.
        """
        self._con.commit()

    def __getstate__(self):
        """
        Return the state of this object as a dictionary (used for pickling).
        """
        state_dict = self.__dict__.copy()
        del state_dict['_con']
        return state_dict

    def __setstate__(self, state_dict):
        """
        Set the state of this object from a dictionary (used for pickling).
        """
        self.__dict__.update(state_dict)
        self._con = sqlite3.connect(self._filename)

class UndirectedGraphFileView(GraphFileView):
    """
    Undirected graph database file view
    """
    def _encoding_to_graph(self, enc):
        """
        Return the undirected graph corresponding to the given plain text graph6 encoding.
        """
        return undirected_graph_from_encoding(enc)

    def _graph_to_encoding(self, g):
        """
        Return the plain text graph6 encoding of the given undirected graph.
        """
        return undirected_graph_to_encoding(g)

    def _encoding_length(self):
        """
        Return the length of the graph6 encoding of each undirected graph in this database file.
        """
        assert self._num_vertices <= 62
        num_bits = (self._num_vertices*(self._num_vertices-1)) // 2
        encoding_length = 1 + ceil(num_bits / 6.0) # graph6 length
        return encoding_length

class DirectedGraphFileView(GraphFileView):
    """
    Directed graph database file view
    """
    def _encoding_to_graph(self, enc):
        """
        Return the directed graph corresponding to the given plain text digraph6 encoding.
        """
        return directed_graph_from_encoding(enc)

    def _graph_to_encoding(self, g):
        """
        Return the plain text digraph6 encoding of the given directed graph.
        """
        return directed_graph_to_encoding(g)

    def _encoding_length(self):
        """
        Return the length of the digraph6 encoding of each directed graph in this database file.
        """
        assert self._num_vertices <= 62
        num_bits = self._num_vertices * self._num_vertices
        encoding_length = len('&') + 1 + ceil(num_bits / 6.0) # digraph6 length
        return encoding_length

class UndirectedToDirectedGraphFileView:
    """
    Undirected to directed graph database file view
    """
    def __init__(self, filename):
        """
        Initialize this "undirected to directed graph" database file view.

        INPUT:

        - ``filename`` -- a string, the path to an SQLite database file (the file will be created if it does not yet exist)
        """
        self._filename = filename
        self._con = sqlite3.connect(filename)
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(name) FROM sqlite_master WHERE type = "table" AND name = "undirected_to_directed"')
        if result.fetchone()[0] == 0:
            cur.execute('CREATE TABLE undirected_to_directed (undirected_graph_id INTEGER, directed_graph_id INTEGER, coefficient INTEGER)')
            cur.execute('CREATE INDEX index_undirected ON undirected_to_directed(undirected_graph_id)')
            self._con.commit()

    def append(self, row):
        """
        Insert a row into the "undirected to directed graph" database file.
        """
        g_idx, h_idx, c = row
        cur = self._con.cursor()
        cur.execute('INSERT INTO undirected_to_directed (undirected_graph_id, directed_graph_id, coefficient) VALUES (?, ?, ?)', (1 + g_idx, 1 + h_idx, c))

    def commit(self):
        """
        Commit any changes to this database file.
        """
        self._con.commit()

    def __iter__(self):
        """
        Return an iterator over the rows in the "undirected to directed graph" database file.
        """
        cur = self._con.cursor()
        for row in cur.execute('SELECT undirected_graph_id, directed_graph_id, coefficient FROM undirected_to_directed'):
            yield (row[0] - 1, row[1] - 1, row[2])

    def __len__(self):
        """
        Return the number of rows in the "undirected to directed graph" database file.
        """
        cur = self._con.cursor()
        result = cur.execute('SELECT COUNT(*) FROM undirected_to_directed').fetchone()
        return result[0]

    def undirected_to_directed_coeffs(self, undirected_graph_idx):
        """
        Return an iterator over the ``(directed_graph_idx, coefficient)`` tuples related to the undirected graph with the given index.
        """
        cur = self._con.cursor()
        for row in cur.execute('SELECT directed_graph_id, coefficient FROM undirected_to_directed WHERE undirected_graph_id = ?', (1 + undirected_graph_idx,)):
            yield (row[0] - 1, row[1])

    def __getstate__(self):
        """
        Return the state of this object as a dictionary (used for pickling).
        """
        state_dict = self.__dict__.copy()
        del state_dict['_con']
        return state_dict

    def __setstate__(self, state_dict):
        """
        Set the state of this object from a dictionary (used for pickling).
        """
        self.__dict__.update(state_dict)
        self._con = sqlite3.connect(self._filename)

class FormalityGraphFileView(GraphFileView):
    """
    Formality graph database file view
    """
    def __init__(self, filename, num_ground_vertices, num_aerial_vertices, num_edges):
        """
        Initialize this formality graph database file view.

        INPUT:

        - ``filename`` -- a string, the path to an SQLite database file (the file will be created if it does not yet exist)

        - ``num_ground_vertices`` -- a natural number, the number of ground vertices in each graph

        - ``num_aerial_vertices`` -- a natural number, the number of aerial vertices in each graph

        - ``num_edges`` -- a natural number, the number of edges in each graph
        """
        self._num_ground_vertices = num_ground_vertices
        self._num_aerial_vertices = num_aerial_vertices
        super().__init__(filename, num_ground_vertices + num_aerial_vertices, num_edges)

    def _encoding_to_graph(self, enc):
        """
        Return the formality graph corresponding to the given encoding.
        """
        return formality_graph_from_encoding((self._num_ground_vertices, self._num_aerial_vertices, enc))

    def _graph_to_encoding(self, g):
        """
        Return the encoding of the given formality graph.
        """
        return formality_graph_to_encoding(g)[2]

    def _encoding_length(self):
        """
        Return the length of the encoding of each formality graph in this database file.
        """
        assert self._num_vertices <= 62
        num_bits = self._num_vertices * self._num_vertices
        encoding_length = 1 + ceil(num_bits / 6.0) # dig6 length
        return encoding_length
