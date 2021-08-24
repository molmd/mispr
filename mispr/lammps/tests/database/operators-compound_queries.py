from pprint import pprint

from first_connection import Connect

if __name__ == "__main__":
    # Connect to local database
    connection = Connect.get_connection()

    # Access the test database
    db = connection.test

    # Read data with embedded fields and comparison operators
    print("Read data with embedded fields and comparison operators: ")
    cursor = db.inventory.find({"size.h": {"$lt": 15}})
    for inventory in cursor:
        pprint(inventory)

    # Read data with compound queries
    print("\nRead data with compound queries: ")
    cursor = db.inventory.find({"status": "A", "qty": {"$lt": 30}})
    for inventory in cursor:
        pprint(inventory)

    # Compound query using explicit $and operator
    print("\nCompound query using explicit $and operator: ")
    cursor = db.inventory.find({"$and": [{"status": "A"}, {"qty": {"$lt": 30}}]})
    for inventory in cursor:
        pprint(inventory)

    # Retrieving data with more than one compounding clause
    print("\nRetrieving data with more than one compounding clause")
    cursor = db.inventory.find(
        {"status": "A", "$or": [{"qty": {"$lt": 30}}, {"item": {"$regex": "^p"}}]}
    )
    for inventory in cursor:
        pprint(inventory)
