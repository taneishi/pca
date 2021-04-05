document.addEventListener('DOMContentLoaded', function () {
    d3.csv('https://raw.githubusercontent.com/taneishi/pca/master/iris.csv').then(function(data){
        columns = ['Sepal.Length', 'Sepal.Width', 'Petal.Length', 'Petal.Width'];

        X = data.map(function(row){
            return columns.map(function(col){ return Number(row[col]); })
        });
        
        pc = new PCA().pca(X, 2);

        labels = data.map(function(row){ return row['Species']; });
        pc = pc.map(function(row, i){ return row.concat(labels[i]); });

        labels = Array.from(new Set(labels));

        colors = ['#aea', '#aae', '#eaa']

        datasets = labels.map(function(label){
            return {
                label: label, 
                data: pc.filter(function(row){ return row[4] == label; }),
			    borderColor: '#222',
                backgroundColor: colors[labels.indexOf(label)],
                pointRadius: 5,
                };
        });

		data = {
		  labels: labels,
		  datasets: datasets,
		};

        var mychart = new Chart('container', {
            type: 'scatter',
            data: data,
            options: {
                responsive: false,
                scales: {
                    x: {
                        title: { display: true, text: 'PC1', },
                    },
                    y: {
                        title: { display: true, text: 'PC2', },
                    },
                },
                plugins: {
                    legend: {
                        position: 'top',
                    },
                    title: {
                        display: true,
                        text: 'Principal Component Anaylsis'
                    },
                }
            },
        });
        
    })
});
