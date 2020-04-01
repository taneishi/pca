var mypca = new PCA();
document.addEventListener('DOMContentLoaded', function () {
	var myChart = Highcharts.chart('container', {
		chart: {
			type: 'scatter'
		},
		title: {
			text: 'Principal Component Anaylsis'
		},
		xAxis: {
			categories: ['Apples', 'Bananas', 'Oranges']
		},
		yAxis: {
			title: {
				text: 'Fruit eaten'
			}
		},
		series: [{
			name: 'Jane',
			data: [1, 0, 4]
		}, {
			name: 'John',
			data: [5, 7, 3]
		}]
	});
});
