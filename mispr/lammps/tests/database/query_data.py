from pprint import pprint

from bson.son import SON

from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test

    # Load more data into MongoDB
    db.inventory.insert_many([
        {"item": "journal",
         "qty": 25,
         "size": SON([("h", 14), ("w", 21), ("uom", "cm")]),
         "status": "A"},
        {"item": "notebook",
         "qty": 50,
         "size": SON([("h", 8.5), ("w", 11), ("uom", "in")]),
         "status": "A"},
        {"item": "paper",
         "qty": 100,
         "size": SON([("h", 8.5), ("w", 11), ("uom", "in")]),
         "status": "D"},
        {"item": "planner",
         "qty": 75,
         "size": SON([("h", 22.85), ("w", 30), ("uom", "cm")]),
         "status": "D"},
        {"item": "postcard",
         "qty": 45,
         "size": SON([("h", 10), ("w", 15.25), ("uom", "cm")]),
         "status": "A"}])

    # Retrieve specific documents in a collection
    print("Retrieve specific documents in a collection: ")
    cursor = db.inventory.find({"status": "D"})
    # Iterate over the results
    for inventory in cursor:
        pprint(inventory)

    # Query Data using Embedded Documents as Criteria
    print("\nQuery Data using Embedded Documents as Criteria: ")
    cursor = db.inventory.find({"size": SON([("h", 14), ("w", 21), ("uom", "cm")])})
    # Iterate over the results
    for inventory in cursor:
        pprint(inventory)

    # Query data using embedded documents with dot notation
    print("\nQuery data using embedded documents with dot notation: ")
    cursor = db.inventory.find({"size.uom": "in"})
    for inventory in cursor:
        pprint(inventory)
