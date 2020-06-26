new Vue({
    el: '#app2',
    template: `
        <div class="cascading-dropdown">
            <div class="dropdown">
                <span>Continent:</span>
                <select v-model="selectedContinent">
                    <option value="">Select a Continent</option>
                    <option v-for="(country_obj, country) in places" :value="country">{{country}}</option>
                </select>
            </div>
            <div class="dropdown">
                <span>Country:</span>
                <select :disabled="countries.length == 0" v-model="selectedCountry">
                    <option value="">Select a Country</option>
                    <option v-for="(city_obj, city) in countries">{{city}}</option>
                </select>
            </div>
            <div class="dropdown">
                <span>City:</span>
                <select :disabled="cities.length == 0" v-model="selectedCity">
                    <option value="">Select a City</option>
                    <option v-for="city in cities">{{city}}</option>
                </select>
            </div>
        </div>
    `,
    data: function() {
        return {
            places: {
                "Asia": {
                    "China": ["Beijing", "Shanghai", "Guangzhou", "Tianjin"],
                    "India": ["New Delhi", "Mumbai", "Bangalore", "Chennai"],
                    "Japan": ["Tokyo", "Kyoto", "Nagoya", "Hiroshima"],
                    "Singapore": ["Singapore"],
                    "Malaysia": ["Kuala Lumpur", "Johor Bahru", "George Town", "Kota Kinabalu"]
                },
                "Europe": {
                    "Germany": ["Berlin", "Hamburg", "Munich", "Cologne", "Frankfurt", "Stuttgart"],
                    "UK": ["London", "Birmingham", "Liverpool", "Bristol"],
                    "France": ["Paris", "Marseille", "Bordeaux", "Toulouse"]
                }
            },
            countries: [],
            cities: [],
            selectedContinent: "",
            selectedCountry: "",
            selectedCity: ""
        }
    },
    watch: {
        selectedContinent: function() {
            // Clear previously selected values
            this.countries = [];
            this.cities = [];
            this.selectedCountry = "";
            this.selectedCity = "";
            // Populate list of countries in the second dropdown
            if (this.selectedContinent.length > 0) {
                this.countries = this.places[this.selectedContinent]
            }
        },
        selectedCountry: function() {
            // Clear previously selected values
            this.cities = [];
            this.selectedCity = "";
            // Now we have a continent and country. Populate list of cities in the third dropdown
            if (this.selectedCountry.length > 0) {
                this.cities = this.places[this.selectedContinent][this.selectedCountry]
            }
        }
    }
})