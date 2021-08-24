from pprint import pprint

from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test

    # Retrieve all documents in the inventory collection
    cursor = db.inventory.find({})
    # Iterate over the results
    for inventory in cursor:
        pprint(inventory)
